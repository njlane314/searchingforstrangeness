#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "ImagePipeline/ImageCentering.h"
#include "ImagePipeline/ImageProduction.h"
#include "ImagePipeline/ImageWindowGeometry.h"
#include "ImagePipeline/SemanticClassifier.h"
#include "Products/ImageFeatures.h"

#include <TVector3.h>
#include <cetlib_except/exception.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <ios>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using image::ImageFeatures;

namespace {

std::string trim(std::string value) {
    auto const not_space = [](unsigned char c) {
        return !std::isspace(c);
    };
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(), not_space));
    value.erase(
        std::find_if(value.rbegin(), value.rend(), not_space).base(),
        value.end());
    return value;
}

} // namespace

class ImageProducer : public art::EDProducer {
  public:
    explicit ImageProducer(fhicl::ParameterSet const &pset);
    void produce(art::Event &e) override;

  private:
    art::InputTag fPFPproducer;
    art::InputTag fSLCproducer;
    art::InputTag fHITproducer;
    art::InputTag fWIREproducer;
    art::InputTag fMCPproducer;
    art::InputTag fBKTproducer;
    art::InputTag fSPproducer;
    art::InputTag fVTXproducer;

    bool fIsData{false};

    std::string fBadChannelFile;
    std::set<unsigned int> fBadChannels;

    std::unique_ptr<image::SemanticClassifier> fSemantic;

    const geo::GeometryCore *fGeo{nullptr};
    const detinfo::DetectorProperties *fDetp{nullptr};
    std::unique_ptr<image::ImageCentering> fCentering;
    std::unique_ptr<image::ImageWindowGeometry> fWindowGeometry;

    struct NeutrinoReco {
        std::vector<art::Ptr<recob::Hit>> slice_hits;
        std::optional<TVector3> vertex;
    };

    void loadBadChannels(const std::string &filename);
    std::vector<art::Ptr<recob::Hit>> collectEventHits(const art::Event &event) const;
    NeutrinoReco collectNeutrinoReco(const art::Event &event) const;
};

ImageProducer::ImageProducer(fhicl::ParameterSet const &pset)
    : fPFPproducer{pset.get<art::InputTag>("PFPproducer")}
    , fSLCproducer{pset.get<art::InputTag>("SLCproducer")}
    , fHITproducer{pset.get<art::InputTag>("HITproducer")}
    , fWIREproducer{pset.get<art::InputTag>("WIREproducer")}
    , fMCPproducer{pset.get<art::InputTag>("MCPproducer")}
    , fBKTproducer{pset.get<art::InputTag>("BKTproducer")}
    , fSPproducer{pset.get<art::InputTag>("SPproducer")}
    , fVTXproducer{pset.get<art::InputTag>("VTXproducer", fPFPproducer)}
    , fIsData{pset.get<bool>("IsData", false)}
    , fBadChannelFile{pset.get<std::string>("BadChannelFile", "")}
{
    if (!fBadChannelFile.empty())
        loadBadChannels(fBadChannelFile);

    fGeo = art::ServiceHandle<geo::Geometry>()->provider();
    fDetp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();

    auto const *det_prop_data = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->provider();
    auto const *clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->provider();

    double const tick_period = clock_data->TPCClock().TickPeriod();
    double const drift_velocity = det_prop_data->DriftVelocity();
    image::ImageWindowGeometryConfig const window_config{
        512U,
        512U,
        tick_period * drift_velocity * 1.0e1,
        fGeo->WirePitch(geo::kU),
        fGeo->WirePitch(geo::kV),
        fGeo->WirePitch(geo::kW)};
    fCentering = std::make_unique<image::ImageCentering>(
        fHITproducer, fSPproducer);
    fWindowGeometry =
        std::make_unique<image::ImageWindowGeometry>(
            *fGeo, window_config);

    fSemantic = std::make_unique<image::SemanticClassifier>(fMCPproducer);

    produces<std::vector<ImageFeatures>>("CroppedWindow");
    produces<std::vector<ImageFeatures>>("FullWindow");
}

void ImageProducer::loadBadChannels(const std::string &filename) {
    fBadChannels.clear();
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw art::Exception(art::errors::Configuration) << "Cannot open bad channel file: " << filename;
    }
    auto parse_channel = [&](std::string const &token,
                             std::size_t line_number) -> unsigned int {
        if (token.empty() || token.front() == '-' || token.front() == '+') {
            throw art::Exception(art::errors::Configuration)
                << "Invalid channel token '" << token << "' in "
                << filename << ":" << line_number;
        }

        std::size_t parsed = 0U;
        unsigned long long value = 0U;
        try {
            value = std::stoull(token, &parsed, 10);
        } catch (std::exception const &) {
            throw art::Exception(art::errors::Configuration)
                << "Invalid channel token '" << token << "' in "
                << filename << ":" << line_number;
        }
        if (parsed != token.size() ||
            value > std::numeric_limits<unsigned int>::max()) {
            throw art::Exception(art::errors::Configuration)
                << "Channel token '" << token << "' is out of range in "
                << filename << ":" << line_number;
        }
        return static_cast<unsigned int>(value);
    };

    std::string line;
    std::size_t line_number = 0U;
    while (std::getline(in, line)) {
        ++line_number;
        auto const comment = line.find('#');
        if (comment != std::string::npos)
            line.erase(comment);
        line = trim(std::move(line));
        if (line.empty())
            continue;

        std::stringstream ss(line);
        std::string first_token;
        std::string second_token;
        std::string extra_token;
        ss >> first_token;
        ss >> second_token;
        ss >> extra_token;
        if (!extra_token.empty()) {
            throw art::Exception(art::errors::Configuration)
                << "Expected one channel or an inclusive channel range in "
                << filename << ":" << line_number << ", got '" << line
                << "'.";
        }

        unsigned int const first =
            parse_channel(first_token, line_number);
        unsigned int const second = second_token.empty()
                                        ? first
                                        : parse_channel(second_token,
                                                        line_number);
        if (second < first) {
            throw art::Exception(art::errors::Configuration)
                << "Descending bad-channel range " << first << " "
                << second << " in " << filename << ":" << line_number;
        }

        for (std::uint64_t channel = first; channel <= second; ++channel) {
            fBadChannels.insert(static_cast<unsigned int>(channel));
        }
    }
}

std::vector<art::Ptr<recob::Hit>> ImageProducer::collectEventHits(const art::Event &event) const {
    auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);

    std::vector<art::Ptr<recob::Hit>> hits;
    hits.reserve(hit_handle->size());

    for (std::size_t i = 0; i < hit_handle->size(); ++i) {
        hits.emplace_back(hit_handle, i);
    }

    return hits;
}

ImageProducer::NeutrinoReco
ImageProducer::collectNeutrinoReco(const art::Event &event) const {
    NeutrinoReco result;
    auto pfp_handle = event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);

    std::optional<size_t> nu_index;
    for (size_t i = 0; i < pfp_handle->size(); ++i) {
        auto const &p = pfp_handle->at(i);
        if (!p.IsPrimary())
            continue;
        int pdg = std::abs(p.PdgCode());
        if (pdg == 12 || pdg == 14 || pdg == 16) {
            nu_index = i;
            break;
        }
    }
    if (!nu_index)
        return result;

    art::FindManyP<recob::Slice> pfp_to_slice(pfp_handle, event, fPFPproducer);
    auto slice_handle = event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    art::FindManyP<recob::Hit> slice_to_hits(slice_handle, event, fSLCproducer);

    auto const slices = pfp_to_slice.at(*nu_index);
    if (!slices.empty() && slices.front()) {
        auto const hits = slice_to_hits.at(slices.front().key());
        result.slice_hits.insert(result.slice_hits.end(),
                                 hits.begin(), hits.end());
    }

    try {
        art::FindManyP<recob::Vertex> pfp_to_vtx(pfp_handle, event, fVTXproducer);
        auto const vtxs = pfp_to_vtx.at(*nu_index);
        if (!vtxs.empty() && vtxs.front()) {
            double xyz[3] = {0., 0., 0.};
            vtxs.front()->XYZ(xyz);
            TVector3 vertex(xyz[0], xyz[1], xyz[2]);
            if (std::isfinite(vertex.X()) &&
                std::isfinite(vertex.Y()) &&
                std::isfinite(vertex.Z())) {
                result.vertex = vertex;
            }
        }
    } catch (cet::exception const &ex) {
        mf::LogDebug("ImageProducer") << "PFParticle->Vertex association unavailable for VTXproducer='"
                                      << fVTXproducer.label() << "': " << ex;
    }

    return result;
}

void ImageProducer::produce(art::Event &event) {
    mf::LogDebug("ImageProducer")
        << "Starting image production for event " << event.id();

    auto neutrino = collectNeutrinoReco(event);
    auto event_hits = collectEventHits(event);

    auto remove_unusable_hits =
        [&](std::vector<art::Ptr<recob::Hit>> &hits) {
            hits.erase(
                std::remove_if(
                    hits.begin(), hits.end(),
                    [&](auto const &hit) {
                        return !hit ||
                               fBadChannels.count(hit->Channel()) != 0U;
                    }),
                hits.end());
        };
    remove_unusable_hits(neutrino.slice_hits);
    remove_unusable_hits(event_hits);

    image::PixelImageOptions opts;
    opts.producers = {fWIREproducer, fHITproducer, fMCPproducer, fBKTproducer};
    opts.semantic = fIsData ? nullptr : fSemantic.get();

    mf::LogDebug("ImageProducer")
        << "IsData=" << std::boolalpha << fIsData << ", semantic images "
        << (fIsData ? "DISABLED" : "ENABLED")
        << ", raw sparse features = [adc, nu-slice]";

    image::ImageCenter const center = fCentering->compute(
        event, neutrino.slice_hits, neutrino.vertex,
        fWindowGeometry->trimmingRadius());
    char const *seed_name = "origin fallback";
    if (center.seed == image::ImageCenterSeed::Vertex)
        seed_name = "vertex";
    else if (center.seed ==
             image::ImageCenterSeed::WeightedCentroid)
        seed_name = "weighted spacepoint centroid";

    mf::LogDebug("ImageProducer")
        << "Cropped-window center = (" << center.position.X() << ", "
        << center.position.Y() << ", " << center.position.Z()
        << "), seed=" << seed_name;

    image::RasterizedImageWindow cropped_window(
        fWindowGeometry->croppedWindowProperties(
            center.position));
    image::RasterizedImageWindow full_window(
        fWindowGeometry->fullWindowProperties());
    image::ImageEventContext context(
        event, event_hits, neutrino.slice_hits, opts);

    std::vector<image::RasterizedImageWindow *> windows{&full_window};
    if (!neutrino.slice_hits.empty())
        windows.insert(windows.begin(), &cropped_window);
    else
        mf::LogDebug("ImageProducer")
            << "No neutrino slice hits found; storing three metadata-bearing "
               "empty CroppedWindow planes for event "
            << event.id();

    image::ImageRasterizer(*fGeo).rasterize(context, windows, fDetp);

    auto cropped = image::SparseImagePacker::pack(
        cropped_window, neutrino.vertex, !fIsData);
    auto full = image::SparseImagePacker::pack(
        full_window, neutrino.vertex, !fIsData);

    auto const cropped_validation =
        image::validateImagePlaneTriplet(cropped);
    if (!cropped_validation) {
        throw cet::exception("ImageProducer")
            << "Internal CroppedWindow contract violation for event "
            << event.id() << ": " << cropped_validation.error;
    }
    auto const full_validation =
        image::validateImagePlaneTriplet(full);
    if (!full_validation) {
        throw cet::exception("ImageProducer")
            << "Internal FullWindow contract violation for event "
            << event.id() << ": " << full_validation.error;
    }

    auto out_cropped =
        std::make_unique<std::vector<ImageFeatures>>(std::move(cropped));
    auto out_full =
        std::make_unique<std::vector<ImageFeatures>>(std::move(full));
    auto const cropped_size = out_cropped->size();
    auto const full_size = out_full->size();
    event.put(std::move(out_cropped), "CroppedWindow");
    event.put(std::move(out_full), "FullWindow");

    mf::LogInfo("ImageProducer")
        << "Stored " << cropped_size << " CroppedWindow and "
        << full_size << " FullWindow ImageFeatures for event "
        << event.id();
}

DEFINE_ART_MODULE(ImageProducer)
