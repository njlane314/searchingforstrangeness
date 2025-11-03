#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Provenance/ProductID.h"
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

#include "sbndcode/BlipReco/Alg/BlipRecoAlg.h"

using image::Image;
using image::ImageCentering;
using image::ImageProperties;
using image::ImageProduct;
using image::SemanticClassifier;

class ImageProducer : public art::EDProducer {
  public:
    explicit ImageProducer(fhicl::ParameterSet const &p);
    void produce(art::Event &e) override;

  private:
    art::InputTag fPFPproducer;
    art::InputTag fSLCproducer;
    art::InputTag fHITproducer;
    art::InputTag fWIREproducer;
    art::InputTag fMCPproducer;
    art::InputTag fBKTproducer;
    bool fIsData{false};

    art::InputTag fT0producer;
    std::unique_ptr<calo::CalorimetryAlg> fCalo;

    std::unique_ptr<blip::BlipRecoAlg> fBlipAlg;
    fhicl::ParameterSet                fBlipAlgPSet;
    art::InputTag                      fBlipAlgHitProducer;

    int fImgW{512};
    int fImgH{512};
    float fADCThresh{4.0f};
    float fCentroidRadiusCm{0.f};

    std::set<unsigned int> fBadChannels;
    std::string fBadChannelFile;

    const geo::GeometryCore *fGeo{nullptr};
    const detinfo::DetectorProperties *fDetp{nullptr};
    double fDriftStepCm{0.0};
    double fPitchU{0.0};
    double fPitchV{0.0};
    double fPitchW{0.0};

    std::unique_ptr<SemanticClassifier> fSemantic;

    void loadBadChannels(const std::string &filename);
    static std::vector<art::Ptr<recob::Hit>>
    collectAllHits(const art::Event &event, art::InputTag const &hitTag);
    std::vector<art::Ptr<recob::Hit>>
    collectNeutrinoSliceHits(const art::Event &event) const;

    std::optional<std::pair<art::ProductID, std::vector<uint8_t>>>
    buildBlipMask(art::Event &event);
};

ImageProducer::ImageProducer(fhicl::ParameterSet const &p) {
    fPFPproducer = p.get<art::InputTag>("PFPproducer");
    fSLCproducer = p.get<art::InputTag>("SLCproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fWIREproducer = p.get<art::InputTag>("WIREproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fBKTproducer = p.get<art::InputTag>("BKTproducer");
    fT0producer = p.get<art::InputTag>("T0producer", art::InputTag{});
    fIsData = p.get<bool>("IsData", false);
    if (p.has_key("BlipAlg")) {
        fBlipAlgPSet = p.get<fhicl::ParameterSet>("BlipAlg");
        fBlipAlg = std::make_unique<blip::BlipRecoAlg>(fBlipAlgPSet);
        fBlipAlgHitProducer = fBlipAlgPSet.get<art::InputTag>("HitProducer", fHITproducer);
    }

    if (p.has_key("CalorimetryAlg")) {
        auto cset = p.get<fhicl::ParameterSet>("CalorimetryAlg");
        fCalo = std::make_unique<calo::CalorimetryAlg>(cset);
    }

    fBadChannelFile = p.get<std::string>("BadChannelFile", "");
    if (!fBadChannelFile.empty())
        loadBadChannels(fBadChannelFile);

    fImgW = p.get<int>("ImageWidth", 512);
    fImgH = p.get<int>("ImageHeight", 512);
    fADCThresh = p.get<float>("ADCImageThreshold", 4.0);
    const float pixel_size_cm = 0.3f;
    fCentroidRadiusCm = 0.5f * std::min(fImgW, fImgH) * pixel_size_cm;

    fGeo = art::ServiceHandle<geo::Geometry>()->provider();
    fDetp =
        art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    auto clock =
        art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
    double tick_period = clock->TPCClock().TickPeriod();
    double drift_vel = fDetp->DriftVelocity();
    fDriftStepCm = tick_period * drift_vel * 1.0e1;
    fPitchU = fGeo->WirePitch(geo::kU);
    fPitchV = fGeo->WirePitch(geo::kV);
    fPitchW = fGeo->WirePitch(geo::kW);

    fSemantic = std::make_unique<SemanticClassifier>(fMCPproducer);

    produces<std::vector<ImageProduct>>("primary_slice");
    produces<std::vector<ImageProduct>>("all_hits");
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
                              art::InputTag const &hitTag) {
    std::vector<art::Ptr<recob::Hit>> out;
    auto h = event.getValidHandle<std::vector<recob::Hit>>(hitTag);
    out.reserve(h->size());
    for (size_t i = 0; i < h->size(); ++i)
        out.emplace_back(h, i);
    return out;
}

std::vector<art::Ptr<recob::Hit>>
ImageProducer::collectNeutrinoSliceHits(const art::Event &event) const {
    auto pfpHandle =
        event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::FindManyP<recob::Slice> pfpToSlice(pfpHandle, event, fPFPproducer);

    std::optional<size_t> pfpIndex;
    for (size_t i = 0; i < pfpHandle->size(); ++i) {
        auto const &p = pfpHandle->at(i);
        if (!p.IsPrimary())
            continue;
        int pdg = std::abs(p.PdgCode());
        if (pdg == 12 || pdg == 14 || pdg == 16) {
            pfpIndex = i;
            break;
        }
    }
    if (!pfpIndex) {
        for (size_t i = 0; i < pfpHandle->size(); ++i) {
            if (pfpHandle->at(i).IsPrimary()) {
                pfpIndex = i;
                break;
            }
        }
    }

    std::vector<art::Ptr<recob::Hit>> result;
    auto sliceHandle =
        event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    art::FindManyP<recob::Hit> sliceToHits(sliceHandle, event, fSLCproducer);

    if (pfpIndex) {
        auto slices = pfpToSlice.at(*pfpIndex);
        if (!slices.empty()) {
            auto const &sl = slices.front();
            auto hits = sliceToHits.at(sl.key());
            result.insert(result.end(), hits.begin(), hits.end());
            return result;
        }
    }

    if (!sliceHandle->empty()) {
        auto hits = sliceToHits.at(0);
        result.insert(result.end(), hits.begin(), hits.end());
    }
    return result;
}

std::optional<std::pair<art::ProductID, std::vector<uint8_t>>>
ImageProducer::buildBlipMask(art::Event &event) {
    if (!fBlipAlg)
        return std::nullopt;

    fBlipAlg->RunBlipReco(event);

    auto hitH = event.getHandle<std::vector<recob::Hit>>(fBlipAlgHitProducer);
    if (!hitH)
        return std::nullopt;

    const auto &hi = fBlipAlg->hitinfo;
    const size_t N = std::min(hi.size(), hitH->size());
    if (N == 0)
        return std::nullopt;

    std::vector<uint8_t> isBlip(N, 0);
    for (size_t i = 0; i < N; ++i) {
        if (hi[i].blipid >= 0)
            isBlip[i] = 1;
    }
    return std::make_pair(hitH.id(), std::move(isBlip));
}

void ImageProducer::produce(art::Event &event) {
    auto all_hits = collectAllHits(event, fHITproducer);
    auto neutrino_hits = collectNeutrinoSliceHits(event);

    auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const detPropDAT =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);

    double T0_ticks = 0.0;
    if (fCalo && !fT0producer.label().empty()) {
        auto pfpHandle =
            event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
        art::FindManyP<recob::Slice> pfpToSlice(pfpHandle, event, fPFPproducer);

        std::optional<size_t> pfpIndex;
        for (size_t i = 0; i < pfpHandle->size(); ++i) {
            auto const &p = pfpHandle->at(i);
            if (!p.IsPrimary()) continue;
            int pdg = std::abs(p.PdgCode());
            if (pdg == 12 || pdg == 14 || pdg == 16) { pfpIndex = i; break; }
        }
        if (!pfpIndex) {
            for (size_t i = 0; i < pfpHandle->size(); ++i)
                if (pfpHandle->at(i).IsPrimary()) { pfpIndex = i; break; }
        }
        if (pfpIndex) {
            auto sliceHandle =
                event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
            art::FindManyP<anab::T0> slcToT0(sliceHandle, event, fT0producer);
            auto slices = pfpToSlice.at(*pfpIndex);
            if (!slices.empty()) {
                auto t0s = slcToT0.at(slices.front().key());
                if (!t0s.empty()) {
                    const double T0_us = t0s.front()->Time();
                    T0_ticks = T0_us / clockData.TPCClock().TickPeriod();
                }
            }
        }
    }

    if (fBlipAlg) {
        if (auto mask = buildBlipMask(event)) {
            const art::ProductID pid = mask->first;
            const auto &isBlip = mask->second;
            auto apply_mask = [&](std::vector<art::Ptr<recob::Hit>> &v) {
                v.erase(std::remove_if(v.begin(), v.end(),
                                       [&](auto const &h) {
                                           const auto same = (h.id() == pid);
                                           const auto idx = static_cast<size_t>(h.key());
                                           return same && idx < isBlip.size() && isBlip[idx];
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
        auto pfpHandle =
            event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
        art::FindManyP<recob::Vertex> pfpToVtx(pfpHandle, event, fPFPproducer);
        for (size_t i = 0; i < pfpHandle->size(); ++i) {
            auto const &p = pfpHandle->at(i);
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

    const double R_U = 0.5 * std::min(fImgH * fDriftStepCm, fImgW * fPitchU);
    const double R_V = 0.5 * std::min(fImgH * fDriftStepCm, fImgW * fPitchV);
    const double R_W = 0.5 * std::min(fImgH * fDriftStepCm, fImgW * fPitchW);

    auto cU = ImageCentering::centroidWithinRadius(
        event, common::TPC_VIEW_U, neutrino_hits, R_U, fBadChannels,
        vtxU.Z(), vtxU.X());
    auto cV = ImageCentering::centroidWithinRadius(
        event, common::TPC_VIEW_V, neutrino_hits, R_V, fBadChannels,
        vtxV.Z(), vtxV.X());
    auto cW = ImageCentering::centroidWithinRadius(
        event, common::TPC_VIEW_W, neutrino_hits, R_W, fBadChannels,
        vtxW.Z(), vtxW.X());

    auto fused = image::fuse_and_project(cU.second, cU.first,
                                         cV.second, cV.first,
                                         cW.second, cW.first);

    std::vector<ImageProperties> props;
    props.emplace_back(fused.wU_star, fused.x_star,
                       fImgW, fImgH, fDriftStepCm, fPitchU,
                       geo::kU);
    props.emplace_back(fused.wV_star, fused.x_star,
                       fImgW, fImgH, fDriftStepCm, fPitchV, geo::kV);
    props.emplace_back(fused.wW_star, fused.x_star,
                       fImgW, fImgH, fDriftStepCm, fPitchW, geo::kW);

    std::vector<Image<float>> det_slice;
    std::vector<Image<int>> sem_slice;
    std::vector<Image<float>> det_event;
    std::vector<Image<int>> sem_event;

    image::PixelImageOptions opts;
    opts.is_data       = fIsData;
    opts.producers     = {fWIREproducer, fHITproducer, fMCPproducer, fBKTproducer};
    opts.adc_threshold = fADCThresh;
    opts.bad_channels  = &fBadChannels;
    opts.semantic      = fSemantic.get();

    image::ImageProduction builder(*fGeo, opts);

    std::optional<image::CalibrationContext> cal;
    if (fCalo) {
        cal = image::CalibrationContext{
            fCalo.get(), &clockData, &detPropDAT, T0_ticks
        };
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

    event.put(std::move(out_slice));
    event.put(std::move(out_event));
}

DEFINE_ART_MODULE(ImageProducer)
