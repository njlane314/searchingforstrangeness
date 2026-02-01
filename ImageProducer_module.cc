#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"

#include "Common/PandoraUtilities.h"
#include "Imaging/Image.h"
#include "Imaging/ImageCentering.h"
#include "Imaging/ImageProduction.h"
#include "Imaging/SemanticClassifier.h"
#include "Products/ImageProducts.h"

#include <TVector3.h>
#include <cetlib_except/exception.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <ios>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using image::Image;
using image::ImageProduct;
using image::ImageProperties;

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

    int fImgW{512};
    int fImgH{512};
    const geo::GeometryCore *fGeo{nullptr};
    const detinfo::DetectorProperties *fDetp{nullptr};
    double fDriftStep{0.0};
    double fPitchU{0.0};
    double fPitchV{0.0};
    double fPitchW{0.0};

    void loadBadChannels(const std::string &filename);
    std::vector<art::Ptr<recob::Hit>> collectEventHits(const art::Event &event) const;
    std::vector<art::Ptr<recob::Hit>> collectNeutrinoSliceHits(const art::Event &event) const;
    std::optional<TVector3> findNeutrinoVertex(const art::Event &event) const;
    TVector3 computeImageCenter(const art::Event &event, const std::vector<art::Ptr<recob::Hit>> &neutrino_hits,
                                std::optional<TVector3> const &vertex, bool &center_defaulted) const;
};

ImageProducer::ImageProducer(fhicl::ParameterSet const &pset) {
    fPFPproducer = pset.get<art::InputTag>("PFPproducer");
    fSLCproducer = pset.get<art::InputTag>("SLCproducer");
    fHITproducer = pset.get<art::InputTag>("HITproducer");
    fWIREproducer = pset.get<art::InputTag>("WIREproducer");
    fMCPproducer = pset.get<art::InputTag>("MCPproducer");
    fBKTproducer = pset.get<art::InputTag>("BKTproducer");
    fSPproducer = pset.get<art::InputTag>("SPproducer");
    fVTXproducer = pset.get<art::InputTag>("VTXproducer", fPFPproducer);

    fIsData = pset.get<bool>("IsData", false);

    fBadChannelFile = pset.get<std::string>("BadChannelFile", "");
    if (!fBadChannelFile.empty())
        loadBadChannels(fBadChannelFile);

    fGeo = art::ServiceHandle<geo::Geometry>()->provider();
    fDetp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();

    auto const *det_prop_data = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->provider();
    auto const *clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->provider();

    double const tick_period = clock_data->TPCClock().TickPeriod();
    double const drift_velocity = det_prop_data->DriftVelocity();
    fDriftStep = tick_period * drift_velocity * 1.0e1;
    fPitchU = fGeo->WirePitch(geo::kU);
    fPitchV = fGeo->WirePitch(geo::kV);
    fPitchW = fGeo->WirePitch(geo::kW);

    fSemantic = std::make_unique<image::SemanticClassifier>(fMCPproducer);

    produces<std::vector<ImageProduct>>("NuSlice");
}

void ImageProducer::loadBadChannels(const std::string &filename) {
    fBadChannels.clear();
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw art::Exception(art::errors::Configuration) << "Cannot open bad channel file: " << filename;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::stringstream ss(line);
        unsigned first = 0;
        unsigned second = 0;
        ss >> first;
        if (ss >> second) {
            for (unsigned ch = first; ch <= second; ++ch)
                fBadChannels.insert(ch);
        } else {
            fBadChannels.insert(first);
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

std::vector<art::Ptr<recob::Hit>> ImageProducer::collectNeutrinoSliceHits(const art::Event &event) const {
    auto pfp_handle = event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::FindManyP<recob::Slice> pfp_to_slice(pfp_handle, event, fPFPproducer);

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

    std::vector<art::Ptr<recob::Hit>> out;
    auto slice_handle = event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    art::FindManyP<recob::Hit> slice_to_hits(slice_handle, event, fSLCproducer);

    if (nu_index) {
        auto slices = pfp_to_slice.at(*nu_index);
        if (!slices.empty()) {
            auto const &sl = slices.front();
            auto hits = slice_to_hits.at(sl.key());
            out.insert(out.end(), hits.begin(), hits.end());
            return out;
        }
    }

    return out;
}

std::optional<TVector3> ImageProducer::findNeutrinoVertex(const art::Event &event) const {
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
        return std::nullopt;

    try {
        art::FindManyP<recob::Vertex> pfp_to_vtx(pfp_handle, event, fVTXproducer);
        auto const vtxs = pfp_to_vtx.at(*nu_index);
        if (vtxs.empty() || !vtxs.front())
            return std::nullopt;

        double xyz[3] = {0., 0., 0.};
        vtxs.front()->XYZ(xyz);
        TVector3 v(xyz[0], xyz[1], xyz[2]);

        if (!std::isfinite(v.X()) || !std::isfinite(v.Y()) || !std::isfinite(v.Z()))
            return std::nullopt;

        return v;
    } catch (cet::exception const &ex) {
        mf::LogDebug("ImageProducer") << "PFParticle->Vertex association unavailable for VTXproducer='"
                                      << fVTXproducer.label() << "': " << ex;
    }

    return std::nullopt;
}

TVector3 ImageProducer::computeImageCenter(const art::Event &event,
                                           const std::vector<art::Ptr<recob::Hit>> &neutrino_hits,
                                           std::optional<TVector3> const &vertex,
                                           bool &center_defaulted) const {
    std::map<art::Ptr<recob::SpacePoint>, double> sp_charge;

    auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    art::FindManyP<recob::SpacePoint> hit_to_sp(hit_handle, event, fSPproducer);

    for (auto const &h : neutrino_hits) {
        if (!h)
            continue;

        double q = std::max(0.f, h->Integral());
        if (!(q > 0.0))
            continue;

        auto const &sps_for_hit = hit_to_sp.at(h.key());
        if (sps_for_hit.empty())
            continue;

        double q_each = q / static_cast<double>(sps_for_hit.size());

        for (auto const &sp : sps_for_hit) {
            if (!sp)
                continue;
            sp_charge[sp] += q_each;
        }
    }

    std::vector<art::Ptr<recob::SpacePoint>> spacepoints;
    std::vector<double> weights;
    spacepoints.reserve(sp_charge.size());
    weights.reserve(sp_charge.size());

    for (auto const &kv : sp_charge) {
        spacepoints.push_back(kv.first);
        weights.push_back(kv.second);
    }

    TVector3 center_world(0., 0., 0.);
    TVector3 center_reference(0., 0., 0.);
    if (vertex) {
        center_reference = *vertex;
    } else {
        center_defaulted = true;
    }

    double Rc = 0.5 * std::min(fImgH * fDriftStep, fImgW * fPitchW);
    if (Rc > 0.0 && std::isfinite(Rc) && !spacepoints.empty()) {
        center_world = image::trimmedCentroid(spacepoints, weights, center_reference, Rc);
    } else {
        center_world = center_reference;
    }

    if (!std::isfinite(center_world.X()) || !std::isfinite(center_world.Y()) || !std::isfinite(center_world.Z())) {
        center_world.SetXYZ(0., 0., 0.);
        center_defaulted = true;
    }

    return center_world;
}

void ImageProducer::produce(art::Event &event) {
    mf::LogDebug("ImageProducer") << "Starting image production for event " << event.id();

    auto neutrino_hits = collectNeutrinoSliceHits(event);
    auto event_hits = collectEventHits(event);

    if (!fBadChannels.empty()) {
        auto remove_bad_channels = [&](std::vector<art::Ptr<recob::Hit>> &hits) {
            hits.erase(std::remove_if(hits.begin(), hits.end(),
                                      [&](auto const &h) {
                                          if (h.isNull())
                                              return false;
                                          return fBadChannels.count(h->Channel()) > 0;
                                      }),
                       hits.end());
        };
        remove_bad_channels(neutrino_hits);
        remove_bad_channels(event_hits);
    }

    if (neutrino_hits.empty()) {
        mf::LogDebug("ImageProducer") << "No neutrino slice hits found; producing empty NuSlice for event "
                                      << event.id();

        auto out_slice = std::make_unique<std::vector<ImageProduct>>();
        event.put(std::move(out_slice), "NuSlice");

        mf::LogInfo("ImageProducer") << "Stored 0 NuSlice ImageProducts for event " << event.id()
                                     << " (no neutrino slice)";
        return;
    }

    auto vertex = findNeutrinoVertex(event);
    bool center_defaulted = false;
    TVector3 center_world = computeImageCenter(event, neutrino_hits, vertex, center_defaulted);

    std::cout << "[ImageProducer] center_world = (" << center_world.X() << ", " << center_world.Y() << ", "
              << center_world.Z() << ") - " << (center_defaulted ? "NO VERTEX (FALLBACK)" : "VERTEX-SEEDED")
              << std::endl;

    TVector3 cU = common::ProjectToWireView(center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_U);
    TVector3 cV = common::ProjectToWireView(center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_V);
    TVector3 cW = common::ProjectToWireView(center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_W);

    std::vector<ImageProperties> props;
    props.emplace_back(cU.Z(), cU.X(), fImgW, fImgH, fDriftStep, fPitchU, geo::kU);
    props.emplace_back(cV.Z(), cV.X(), fImgW, fImgH, fDriftStep, fPitchV, geo::kV);
    props.emplace_back(cW.Z(), cW.X(), fImgW, fImgH, fDriftStep, fPitchW, geo::kW);

    std::vector<Image<float>> det_slice;
    std::vector<Image<int>> sem_slice;

    image::PixelImageOptions opts;
    opts.producers = {fWIREproducer, fHITproducer, fMCPproducer, fBKTproducer};
    opts.semantic = fIsData ? nullptr : fSemantic.get();

    mf::LogDebug("ImageProducer") << "IsData=" << std::boolalpha << fIsData << ", semantic images "
                                  << (fIsData ? "DISABLED" : "ENABLED");

    image::ImageProduction builder(*fGeo, opts);

    builder.build(event, event_hits, props, det_slice, sem_slice, fDetp);

    mf::LogDebug("ImageProducer") << "Built images: det=" << det_slice.size() << ", sem=" << sem_slice.size();

    for (std::size_t i = 0; i < props.size(); ++i) {
        mf::LogDebug("ImageProducer") << "View " << static_cast<int>(props[i].view()) << " size " << props[i].width()
                                      << "x" << props[i].height();
    }

    std::array<uint32_t, 3> hit_counts{0, 0, 0};
    for (auto const &hit : event_hits) {
        switch (hit->View()) {
        case geo::kU:
            ++hit_counts[0];
            break;
        case geo::kV:
            ++hit_counts[1];
            break;
        case geo::kW:
            ++hit_counts[2];
            break;
        default:
            break;
        }
    }

    auto vertex_pixel = [&](ImageProperties const &p) -> std::pair<int32_t, int32_t> {
        if (!vertex)
            return {-1, -1};
        common::PandoraView pandora_view;
        switch (p.view()) {
        case geo::kU:
            pandora_view = common::TPC_VIEW_U;
            break;
        case geo::kV:
            pandora_view = common::TPC_VIEW_V;
            break;
        case geo::kW:
        case geo::kY:
            pandora_view = common::TPC_VIEW_W;
            break;
        default:
            return {-1, -1};
        }

        auto proj = common::ProjectToWireView(vertex->X(), vertex->Y(), vertex->Z(), pandora_view);
        auto col = p.col(proj.Z());
        auto row = p.row(proj.X());
        if (!col || !row)
            return {-1, -1};
        return {static_cast<int32_t>(*row), static_cast<int32_t>(*col)};
    };

    auto pack_plane = [&](Image<float> const &det, Image<int> const &sem, ImageProperties const &p, bool include_sem,
                          uint32_t hit_count) {
        ImageProduct out;
        out.view = static_cast<int>(p.view());
        out.width = static_cast<uint32_t>(p.width());
        out.height = static_cast<uint32_t>(p.height());
        out.hit_count = hit_count;
        out.origin_x = static_cast<float>(p.origin_x());
        out.origin_y = static_cast<float>(p.origin_y());
        out.pixel_w = static_cast<float>(p.pixel_w());
        out.pixel_h = static_cast<float>(p.pixel_h());
        auto [vertex_row, vertex_col] = vertex_pixel(p);
        out.vertex_row = vertex_row;
        out.vertex_col = vertex_col;
        out.adc = det.data();
        if (include_sem) {
            auto tmp = sem.data();
            out.semantic.assign(tmp.begin(), tmp.end());
        }
        return out;
    };

    auto out_slice = std::make_unique<std::vector<ImageProduct>>();
    out_slice->reserve(3);
    for (std::size_t i = 0; i < 3 && i < det_slice.size(); ++i) {
        out_slice->emplace_back(pack_plane(det_slice[i], sem_slice[i], props[i], !fIsData, hit_counts[i]));
    }

    auto const n_slice = out_slice->size();
    event.put(std::move(out_slice), "NuSlice");

    mf::LogInfo("ImageProducer") << "Stored " << n_slice << " NuSlice ImageProducts for event " << event.id();
}

DEFINE_ART_MODULE(ImageProducer)
