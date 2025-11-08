#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
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
#include "lardataobj/RecoBase/Wire.h"

#include "Common/PandoraUtilities.h"
#include "Imaging/Image.h"
#include "Imaging/ImageCentering.h"
#include "Imaging/ImageProduction.h"
#include "Imaging/SemanticClassifier.h"
#include "Products/ImageProducts.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"

#include <TVector3.h>
#include <algorithm>
#include <cmath>
#include <cetlib_except/exception.h>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using image::Image;
using image::ImageProperties;
using image::ImageProduct;

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
    art::InputTag fT0producer;

    bool fIsData{false};

    std::string fBadChannelFile;
    std::set<unsigned int> fBadChannels;

    std::unique_ptr<calo::CalorimetryAlg> fCalo;
    std::unique_ptr<blip::BlipRecoAlg> fBlipAlg;
    std::unique_ptr<sem::SemanticClassifier> fSemantic;

    int fImgW{512};
    int fImgH{512};
    float fADCThresh{4.0f};
    std::map<geo::View_t, double> fCentroidRadiusCm{};

    const geo::GeometryCore *fGeo{nullptr};
    const detinfo::DetectorProperties *fDetp{nullptr};
    double fDriftStepCm{0.0};
    double fPitchU{0.0};
    double fPitchV{0.0};
    double fPitchW{0.0};

    void loadBadChannels(const std::string &filename);
    static std::vector<art::Ptr<recob::Hit>> collectAllHits(const art::Event &event,
                                                            art::InputTag const &hitProducer);
    std::vector<art::Ptr<recob::Hit>> collectNeutrinoSliceHits(const art::Event &event) const;

    std::map<art::Ptr<recob::Hit>, std::size_t> buildBlipMask(art::Event &event);

    double collectNeutrinoTime(art::Event &event) const;
};

ImageProducer::ImageProducer(fhicl::ParameterSet const &pset) {
    fPFPproducer = pset.get<art::InputTag>("PFPproducer");
    fSLCproducer = pset.get<art::InputTag>("SLCproducer");
    fHITproducer = pset.get<art::InputTag>("HITproducer");
    fWIREproducer = pset.get<art::InputTag>("WIREproducer");
    fMCPproducer = pset.get<art::InputTag>("MCPproducer");
    fBKTproducer = pset.get<art::InputTag>("BKTproducer");
    fT0producer = pset.get<art::InputTag>("T0producer");

    fIsData = pset.get<bool>("IsData", false);

    fBadChannelFile = pset.get<std::string>("BadChannelFile", "");
    if (!fBadChannelFile.empty())
        loadBadChannels(fBadChannelFile);

    fImgW = pset.get<int>("ImageWidth", 512);
    fImgH = pset.get<int>("ImageHeight", 512);
    fADCThresh = pset.get<float>("ADCImageThreshold", 4.0);

    fGeo = art::ServiceHandle<geo::Geometry>()->provider();
    fDetp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();

    auto const* det_prop_data =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->provider();
    auto const* clock_data =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->provider();

    double const tick_period = clock_data->TPCClock().TickPeriod();
    double const drift_velocity = det_prop_data->DriftVelocity();
    fDriftStepCm = tick_period * drift_velocity * 1.0e1;
    fPitchU = fGeo->WirePitch(geo::kU);
    fPitchV = fGeo->WirePitch(geo::kV);
    fPitchW = fGeo->WirePitch(geo::kW);

    auto compute_radius = [&](double pitch) {
        return 0.5 * std::min(fImgH * fDriftStepCm, fImgW * pitch);
    };
    fCentroidRadiusCm[geo::kU] = compute_radius(fPitchU);
    fCentroidRadiusCm[geo::kV] = compute_radius(fPitchV);
    fCentroidRadiusCm[geo::kW] = compute_radius(fPitchW);

    if (pset.has_key("CaloAlg")) {
        auto const calo_alg_config = pset.get<fhicl::ParameterSet>("CaloAlg");
        fCalo = std::make_unique<calo::CalorimetryAlg>(calo_alg_config);
    }
    if (pset.has_key("BlipAlg")) {
        fBlipAlg.reset(
            new blip::BlipRecoAlg(pset.get<fhicl::ParameterSet>("BlipAlg")));
    }

    fSemantic = std::make_unique<sem::SemanticClassifier>(fMCPproducer);

    produces<std::vector<ImageProduct>>("neutrino_slice_images");
    produces<std::vector<ImageProduct>>("interaction_event_images");

}

void ImageProducer::loadBadChannels(const std::string &filename) {
    fBadChannels.clear();
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw art::Exception(art::errors::Configuration)
            << "Cannot open bad channel file: " << filename;
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

std::vector<art::Ptr<recob::Hit>>
ImageProducer::collectAllHits(const art::Event &event,
                              art::InputTag const &hitProducer) {
    std::vector<art::Ptr<recob::Hit>> out;
    auto h = event.getValidHandle<std::vector<recob::Hit>>(hitProducer);
    out.reserve(h->size());
    for (size_t i = 0; i < h->size(); ++i)
        out.emplace_back(h, i);
    return out;
}

std::vector<art::Ptr<recob::Hit>>
ImageProducer::collectNeutrinoSliceHits(const art::Event &event) const {
    auto pfp_handle =
        event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
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
    auto slice_handle =
        event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
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

std::map<art::Ptr<recob::Hit>, std::size_t>
ImageProducer::buildBlipMask(art::Event &event) {
    if (!fBlipAlg) return {};

    fBlipAlg->RunBlipReco(event);

    auto hit_h = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);

    auto const &info = fBlipAlg->hitinfo;
    if (info.size() != hit_h->size()) {
        throw art::Exception(art::errors::LogicError)
            << "BlipRecoAlg hitinfo size (" << info.size()
            << ") does not match hit product size (" << hit_h->size() << ").";
    }

    std::map<art::Ptr<recob::Hit>, std::size_t> blip_hit_to_key;
    for (std::size_t i = 0; i < info.size(); ++i) {
        if (info[i].blipid >= 0) {
            art::Ptr<recob::Hit> ph(hit_h, i);
            blip_hit_to_key.emplace(ph, static_cast<std::size_t>(ph.key()));
        }
    }
    return blip_hit_to_key;
}

double ImageProducer::collectNeutrinoTime(art::Event &event) const
{
    auto const* clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->provider();
    if (fT0producer.label().empty()) return 0.0;

    auto pfp_h = event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);

    std::optional<size_t> nu_index;
    for (size_t i = 0; i < pfp_h->size(); ++i) {
        auto const &p = pfp_h->at(i);
        if (!p.IsPrimary()) continue;
        int const pdg = std::abs(p.PdgCode());
        if (pdg == 12 || pdg == 14 || pdg == 16) {
            nu_index = i;
            break;
        }
    }
    if (nu_index) {
        art::FindManyP<recob::Slice> pfp_to_slice(pfp_h, event, fPFPproducer);
        if (pfp_to_slice.isValid()) {
            auto const slices = pfp_to_slice.at(*nu_index);
            if (!slices.empty()) {
                auto slc_h = event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
                art::FindManyP<anab::T0> slc_to_t0(slc_h, event, fT0producer);
                if (slc_to_t0.isValid()) {
                    auto const &t0s = slc_to_t0.at(slices.front().key());
                    if (!t0s.empty()) {
                        double const T0_ns = t0s.front()->Time();
                        return (T0_ns * 1.0e-3) / clockData->TPCClock().TickPeriod();
                    }
                }
            }
        }
    }
    return 0.0;
}

void ImageProducer::produce(art::Event &event) {
    auto all_hits = collectAllHits(event, fHITproducer);
    auto neutrino_hits = collectNeutrinoSliceHits(event);

    auto const* clock_data =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->provider();
    auto const* det_prop =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->provider();

    double T0_ticks = (fCalo ? collectNeutrinoTime(event) : 0.0);

    if (fBlipAlg) {
        auto blip_hit_to_key = buildBlipMask(event);
        if (!blip_hit_to_key.empty()) {
            auto apply_mask = [&](std::vector<art::Ptr<recob::Hit>> &v) {
                v.erase(std::remove_if(v.begin(), v.end(),
                                       [&](auto const &h) {
                                           if (h.isNull()) return false;
                                           return blip_hit_to_key.find(h) != blip_hit_to_key.end();
                                       }),
                        v.end());
            };
            apply_mask(neutrino_hits);
            apply_mask(all_hits);
        }
    }

    double vtx_x = std::numeric_limits<double>::quiet_NaN();
    double vtx_y = std::numeric_limits<double>::quiet_NaN();
    double vtx_z = std::numeric_limits<double>::quiet_NaN();
    {
        auto pfp_handle =
            event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
        art::FindManyP<recob::Vertex> pfpToVtx(pfp_handle, event, fPFPproducer);
        for (size_t i = 0; i < pfp_handle->size(); ++i) {
            auto const &p = pfp_handle->at(i);
            if (!p.IsPrimary())
                continue;
            auto const &v = pfpToVtx.at(i);
            if (!v.empty()) {
                auto const &pos = v.front()->position();
                vtx_x = pos.X();
                vtx_y = pos.Y();
                vtx_z = pos.Z();
                break;
            }
        }
    }

    TVector3 vtx_world(vtx_x, vtx_y, vtx_z);
    TVector3 vtxU = common::ProjectToWireView(
        vtx_world.X(), vtx_world.Y(), vtx_world.Z(), common::TPC_VIEW_U);
    TVector3 vtxV = common::ProjectToWireView(
        vtx_world.X(), vtx_world.Y(), vtx_world.Z(), common::TPC_VIEW_V);
    TVector3 vtxW = common::ProjectToWireView(
        vtx_world.X(), vtx_world.Y(), vtx_world.Z(), common::TPC_VIEW_W);

    const double R_U = fCentroidRadiusCm.at(geo::kU);
    const double R_V = fCentroidRadiusCm.at(geo::kV);
    const double R_W = fCentroidRadiusCm.at(geo::kW);

    auto cU = image::centroidWithinRadius(
        event, common::TPC_VIEW_U, neutrino_hits, R_U, fBadChannels,
        vtxU.Z(), vtxU.X());
    auto cV = image::centroidWithinRadius(
        event, common::TPC_VIEW_V, neutrino_hits, R_V, fBadChannels,
        vtxV.Z(), vtxV.X());
    auto cW = image::centroidWithinRadius(
        event, common::TPC_VIEW_W, neutrino_hits, R_W, fBadChannels,
        vtxW.Z(), vtxW.X());

    std::vector<ImageProperties> props;
    props.emplace_back(cU.first, cU.second,
                       fImgW, fImgH, fDriftStepCm, 
                       fPitchU, geo::kU);
    props.emplace_back(cV.first, cV.second,
                       fImgW, fImgH, fDriftStepCm, 
                       fPitchV, geo::kV);
    props.emplace_back(cW.first, cW.second,
                       fImgW, fImgH, fDriftStepCm, 
                       fPitchW, geo::kW);

    std::vector<Image<float>> det_slice;
    std::vector<Image<int>> sem_slice;
  
    std::vector<Image<float>> det_event;
    std::vector<Image<int>> sem_event;

    image::PixelImageOptions opts;
    opts.producers = {fWIREproducer, fHITproducer, fMCPproducer, fBKTproducer};
    opts.semantic  = fIsData ? nullptr : fSemantic.get();

    image::ImageProduction builder(*fGeo, opts);

    std::optional<image::CalibrationContext> cal;
    if (fCalo) {
        image::CalibrationContext tmp;
        tmp.calo = fCalo.get();
        tmp.clocks = clock_data;
        tmp.detprop = det_prop;
        tmp.T0_ticks = T0_ticks;
        cal = tmp;
    }

    builder.build(
        event,
        neutrino_hits,
        props,
        det_slice,
        sem_slice,
        fDetp,
        cal);

    builder.build(
        event,
        all_hits,
        props,
        det_event,
        sem_event,
        fDetp,
        cal);

    auto pack_plane = [](Image<float> const &det, Image<int> const &sem,
                         ImageProperties const &p, bool include_sem) {
        ImageProduct out;
        out.view = static_cast<int>(p.view());
        out.width = static_cast<uint32_t>(p.width());
        out.height = static_cast<uint32_t>(p.height());
        out.origin_x = static_cast<float>(p.origin_x());
        out.origin_y = static_cast<float>(p.origin_y());
        out.pixel_w = static_cast<float>(p.pixel_w());
        out.pixel_h = static_cast<float>(p.pixel_h());
        out.adc = det.data();
        if (include_sem) {
            auto tmp = sem.data();
            out.semantic.assign(tmp.begin(), tmp.end());
        }
        return out;
    };

    auto out_slice = std::make_unique<std::vector<ImageProduct>>();
    auto out_event = std::make_unique<std::vector<ImageProduct>>();
    out_slice->reserve(3);
    out_event->reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        out_slice->emplace_back(
            pack_plane(det_slice[i], sem_slice[i], props[i], !fIsData));
        out_event->emplace_back(
            pack_plane(det_event[i], sem_event[i], props[i], !fIsData));
    }

    event.put(std::move(out_slice), "neutrino_slice_images");
    event.put(std::move(out_event), "interaction_event_images");
}

DEFINE_ART_MODULE(ImageProducer)
