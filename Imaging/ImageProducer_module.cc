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
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "Products/ImageProducts.h"
#include "Common/PandoraUtilities.h"
#include "Imaging/Image.h"
#include "Imaging/ImageAlgo.h"
#include "Imaging/SemanticClassifier.h"

#include <TVector3.h>
#include <algorithm>
#include <cmath>
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
using image::ImageAlgo;
using image::ImageProperties;
using image::PlaneImage;
using image::SemanticClassifier;

namespace {
static std::pair<double, double> centroidWithinRadius(const art::Event &event,
                                                       common::PandoraView view,
                                                       const std::vector<art::Ptr<recob::Hit>> &hits,
                                                       double radius,
                                                       const std::set<unsigned int> &bad_channels,
                                                       double vtx_z,
                                                       double vtx_x) {
  double W = 0.0;
  double Zs = 0.0;
  double Xs = 0.0;
  for (auto const &h : hits) {
    if (bad_channels.count(h->Channel())) continue;
    if (common::GetPandoraView(h) != view) continue;
    double q = std::max(0.f, h->Integral());
    if (q <= 0) continue;
    TVector3 p = common::GetPandoraHitPosition(event, h, view);
    double d = std::hypot(p.Z() - vtx_z, p.X() - vtx_x);
    if (d <= radius) {
      W += q;
      Zs += q * p.Z();
      Xs += q * p.X();
    }
  }
  if (W == 0.0) return {vtx_z, vtx_x};
  return {Zs / W, Xs / W};
}
} // namespace

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
  std::unique_ptr<ImageAlgo> fAlgo;

  void loadBadChannels(const std::string &filename);
  static std::vector<art::Ptr<recob::Hit>> collectAllHits(const art::Event &event,
                                                          art::InputTag const &hitTag);
  std::vector<art::Ptr<recob::Hit>> collectNeutrinoSliceHits(const art::Event &event) const;
};

ImageProducer::ImageProducer(fhicl::ParameterSet const &p) {
  fPFPproducer = p.get<art::InputTag>("PFPproducer");
  fSLCproducer = p.get<art::InputTag>("SLCproducer");
  fHITproducer = p.get<art::InputTag>("HITproducer");
  fWIREproducer = p.get<art::InputTag>("WIREproducer");
  fMCPproducer = p.get<art::InputTag>("MCPproducer");
  fBKTproducer = p.get<art::InputTag>("BKTproducer");
  fIsData = p.get<bool>("IsData", false);

  fBadChannelFile = p.get<std::string>("BadChannelFile", "");
  if (!fBadChannelFile.empty()) loadBadChannels(fBadChannelFile);

  fImgW = p.get<int>("ImageWidth", 512);
  fImgH = p.get<int>("ImageHeight", 512);
  fADCThresh = p.get<float>("ADCImageThreshold", 4.0);
  const float pixel_size_cm = 0.3f;
  fCentroidRadiusCm = 0.5f * std::min(fImgW, fImgH) * pixel_size_cm;

  fGeo = art::ServiceHandle<geo::Geometry>()->provider();
  fDetp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
  auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
  double tick_period = clock->TPCClock().TickPeriod();
  double drift_vel = fDetp->DriftVelocity();
  fDriftStepCm = tick_period * drift_vel * 1.0e1;
  fPitchU = fGeo->WirePitch(geo::kU);
  fPitchV = fGeo->WirePitch(geo::kV);
  fPitchW = fGeo->WirePitch(geo::kW);

  fSemantic = std::make_unique<SemanticClassifier>(fMCPproducer);
  fAlgo = std::make_unique<ImageAlgo>(fWIREproducer, fHITproducer, fMCPproducer, fBKTproducer,
                                      fADCThresh, fGeo, fDetp);

  produces<std::vector<PlaneImage>>("slice");
  produces<std::vector<PlaneImage>>("event");
}

void ImageProducer::loadBadChannels(const std::string &filename) {
  fBadChannels.clear();
  std::ifstream in(filename);
  if (!in.is_open()) {
    throw art::Exception(art::errors::Configuration) << "Cannot open bad channel file: " << filename;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line.front() == '#') continue;
    std::stringstream ss(line);
    unsigned first = 0;
    unsigned second = 0;
    ss >> first;
    if (ss >> second) {
      for (unsigned ch = first; ch <= second; ++ch) fBadChannels.insert(ch);
    } else {
      fBadChannels.insert(first);
    }
  }
}


std::vector<art::Ptr<recob::Hit>> ImageProducer::collectAllHits(const art::Event &event,
                                                                art::InputTag const &hitTag) {
  std::vector<art::Ptr<recob::Hit>> out;
  auto h = event.getValidHandle<std::vector<recob::Hit>>(hitTag);
  out.reserve(h->size());
  for (size_t i = 0; i < h->size(); ++i) out.emplace_back(h, i);
  return out;
}

std::vector<art::Ptr<recob::Hit>> ImageProducer::collectNeutrinoSliceHits(const art::Event &event) const {
  auto pfpHandle = event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
  art::FindManyP<recob::Slice> pfpToSlice(pfpHandle, event, fPFPproducer);

  std::optional<size_t> pfpIndex;
  for (size_t i = 0; i < pfpHandle->size(); ++i) {
    auto const &p = pfpHandle->at(i);
    if (!p.IsPrimary()) continue;
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
  auto sliceHandle = event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
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

void ImageProducer::produce(art::Event &event) {
  auto all_hits = collectAllHits(event, fHITproducer);
  auto neutrino_hits = collectNeutrinoSliceHits(event);

  double vtx_x = std::numeric_limits<double>::quiet_NaN();
  double vtx_y = std::numeric_limits<double>::quiet_NaN();
  double vtx_z = std::numeric_limits<double>::quiet_NaN();
  {
    auto pfpHandle = event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::FindManyP<recob::Vertex> pfpToVtx(pfpHandle, event, fPFPproducer);
    for (size_t i = 0; i < pfpHandle->size(); ++i) {
      auto const &p = pfpHandle->at(i);
      if (!p.IsPrimary()) continue;
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
  TVector3 vtxU = common::ProjectToWireView(vtx_world.X(), vtx_world.Y(), vtx_world.Z(), common::TPC_VIEW_U);
  TVector3 vtxV = common::ProjectToWireView(vtx_world.X(), vtx_world.Y(), vtx_world.Z(), common::TPC_VIEW_V);
  TVector3 vtxW = common::ProjectToWireView(vtx_world.X(), vtx_world.Y(), vtx_world.Z(), common::TPC_VIEW_W);

  auto cU = centroidWithinRadius(event, common::TPC_VIEW_U, neutrino_hits, fCentroidRadiusCm, fBadChannels,
                                 vtxU.Z(), vtxU.X());
  auto cV = centroidWithinRadius(event, common::TPC_VIEW_V, neutrino_hits, fCentroidRadiusCm, fBadChannels,
                                 vtxV.Z(), vtxV.X());
  auto cW = centroidWithinRadius(event, common::TPC_VIEW_W, neutrino_hits, fCentroidRadiusCm, fBadChannels,
                                 vtxW.Z(), vtxW.X());

  std::vector<ImageProperties> props;
  props.emplace_back(cU.first, cU.second, fImgW, fImgH, fDriftStepCm, fPitchU, geo::kU);
  props.emplace_back(cV.first, cV.second, fImgW, fImgH, fDriftStepCm, fPitchV, geo::kV);
  props.emplace_back(cW.first, cW.second, fImgW, fImgH, fDriftStepCm, fPitchW, geo::kW);

  std::vector<Image<float>> det_slice;
  std::vector<Image<int>> sem_slice;
  std::vector<Image<float>> det_event;
  std::vector<Image<int>> sem_event;
  fAlgo->produceImages(event, neutrino_hits, props, fIsData, fSemantic.get(), fBadChannels, det_slice, sem_slice);
  fAlgo->produceImages(event, all_hits, props, fIsData, fSemantic.get(), fBadChannels, det_event, sem_event);

  auto pack_plane = [](Image<float> const &det, Image<int> const &sem, ImageProperties const &p, bool include_sem) {
    PlaneImage out;
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

  auto out_slice = std::make_unique<std::vector<PlaneImage>>();
  auto out_event = std::make_unique<std::vector<PlaneImage>>();
  out_slice->reserve(3);
  out_event->reserve(3);
  for (size_t i = 0; i < 3; ++i) {
    out_slice->emplace_back(pack_plane(det_slice[i], sem_slice[i], props[i], !fIsData));
    out_event->emplace_back(pack_plane(det_event[i], sem_event[i], props[i], !fIsData));
  }

  event.put(std::move(out_slice), "slice");
  event.put(std::move(out_event), "event");
}

DEFINE_ART_MODULE(ImageProducer)
