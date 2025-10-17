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
#include "Products/SegmentationProducts.h"
#include "Products/InferencePerf.h"
#include "Common/PandoraUtilities.h"
#include "Imaging/Image.h"
#include "Imaging/ImageAlgo.h"
#include "Imaging/SemanticPixelClassifier.h"

#include <TVector3.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include <limits.h>

using namespace image;

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

  std::string fContainerImage;
  std::string fWrapperScript;
  std::vector<std::string> fWrapperExtra;
  bool fWantClassification{true};
  bool fWantSegmentation{true};

  int fImgW{512};
  int fImgH{512};
  float fADCThresh{4.0f};
  float fCentroidRadiusCm{0.f};

  std::set<unsigned int> fBadChannels;
  std::string fBadChannelFile;

  std::vector<ModelConfig> fModels;
  std::vector<std::string> fActiveModels;
  std::string fAssetsBaseDir;
  std::string fWeightsBaseDir;
  std::string fInferenceWrapper;

  const geo::GeometryCore *fGeo{nullptr};
  const detinfo::DetectorProperties *fDetp{nullptr};
  double fDriftStepCm{0.0};
  double fPitchU{0.0};
  double fPitchV{0.0};
  double fPitchW{0.0};

  std::unique_ptr<SemanticPixelClassifier> fSemantic;
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

  fContainerImage = p.get<std::string>("ContainerImage");
  fWrapperScript = p.get<std::string>("WrapperScript", "/app/infer_bin.py");
  fWrapperExtra = p.get<std::vector<std::string>>("WrapperExtraArgs", {});
  fWantClassification = p.get<bool>("WantClassification", true);
  fWantSegmentation = p.get<bool>("WantSegmentation", true);

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

  fAssetsBaseDir = p.get<std::string>("AssetsBaseDir", "");
  fWeightsBaseDir = p.get<std::string>("WeightsBaseDir", "weights");
  fInferenceWrapper = p.get<std::string>("InferenceWrapper", "scripts/inference_wrapper.sh");

  auto model_psets = p.get<std::vector<fhicl::ParameterSet>>("Models", {});
  for (auto const &ps : model_psets) {
    fModels.push_back({ps.get<std::string>("name"),
                       ps.get<std::string>("weights_file"),
                       ps.get<std::string>("arch")});
  }
  fActiveModels = p.get<std::vector<std::string>>("ActiveModels", {});

  fSemantic = std::make_unique<SemanticPixelClassifier>(fMCPproducer);
  fAlgo = std::make_unique<ImageAlgo>(fWIREproducer, fHITproducer, fMCPproducer, fBKTproducer,
                                      fADCThresh, fWeightsBaseDir, fInferenceWrapper, fAssetsBaseDir,
                                      fModels, fActiveModels, fGeo, fDetp, ".");

  produces<std::vector<PlaneImage>>("slice");
  produces<std::vector<PlaneImage>>("event");
  produces<InferenceScores>();
  produces<std::vector<PlaneSegmentation>>("seg");
  produces<InferencePerfProduct>("perf");
}

void ImageProducerED::loadBadChannels(const std::string &filename) {
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


std::vector<art::Ptr<recob::Hit>> ImageProducerED::collectAllHits(const art::Event &event,
                                                                  art::InputTag const &hitTag) {
  std::vector<art::Ptr<recob::Hit>> out;
  auto h = event.getValidHandle<std::vector<recob::Hit>>(hitTag);
  out.reserve(h->size());
  for (size_t i = 0; i < h->size(); ++i) out.emplace_back(h, i);
  return out;
}

std::vector<art::Ptr<recob::Hit>> ImageProducerED::collectNeutrinoSliceHits(const art::Event &event) const {
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

void ImageProducerED::produce(art::Event &event) {
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

  std::string scratch = ".";
  if (auto *s = std::getenv("_CONDOR_SCRATCH_DIR")) scratch = s;
  char buf[PATH_MAX];
  std::string abs_scratch = realpath(scratch.c_str(), buf) ? std::string(buf) : scratch;

  BinaryInferenceOutput inf_out = fAlgo->runInferenceBinary(det_slice, abs_scratch,
                                                            fWantClassification, fWantSegmentation);

  auto sc = std::make_unique<InferenceScores>();
  for (auto const &kv : inf_out.scores) {
    sc->names.push_back(kv.first);
    sc->scores.push_back(kv.second);
  }

  auto perf_prod = std::make_unique<InferencePerfProduct>();
  perf_prod->per_model.reserve(inf_out.perfs.size());
  for (auto const &kv : inf_out.perfs) {
    ModelPerf mp;
    mp.model = kv.first;
    mp.t_write_req_ms   = static_cast<float>(kv.second.t_write_req_ms);
    mp.t_exec_total_ms  = static_cast<float>(kv.second.t_exec_total_ms);
    mp.t_read_resp_ms   = static_cast<float>(kv.second.t_read_resp_ms);
    mp.t_child_total_ms = static_cast<float>(kv.second.t_child_total_ms);
    mp.t_child_setup_ms = static_cast<float>(kv.second.t_child_setup_ms);
    mp.t_child_infer_ms = static_cast<float>(kv.second.t_child_infer_ms);
    mp.t_child_post_ms  = static_cast<float>(kv.second.t_child_post_ms);
    mp.child_max_rss_mb = static_cast<float>(kv.second.child_max_rss_mb);
    mp.child_cuda_mem_mb= static_cast<float>(kv.second.child_cuda_mem_mb);
    perf_prod->per_model.push_back(std::move(mp));
  }

  std::unique_ptr<std::vector<PlaneSegmentation>> seg_prod;
  if (inf_out.has_segmentation) {
    seg_prod = std::make_unique<std::vector<PlaneSegmentation>>();
    seg_prod->reserve(3);
    auto make_plane = [&](int v, const std::vector<uint8_t> &lbl, const std::vector<float> &conf) {
      PlaneSegmentation p;
      p.view = v;
      p.width = inf_out.segW;
      p.height = inf_out.segH;
      p.nclasses = 0;
      p.labels = lbl;
      p.confidence = conf;
      return p;
    };
    seg_prod->push_back(make_plane(int(geo::kU), inf_out.seg_u, inf_out.seg_conf_u));
    seg_prod->push_back(make_plane(int(geo::kV), inf_out.seg_v, inf_out.seg_conf_v));
    seg_prod->push_back(make_plane(int(geo::kW), inf_out.seg_w, inf_out.seg_conf_w));
  }

  event.put(std::move(out_slice), "slice");
  event.put(std::move(out_event), "event");
  event.put(std::move(sc));
  if (seg_prod) event.put(std::move(seg_prod), "seg");
  event.put(std::move(perf_prod), "perf");
}

DEFINE_ART_MODULE(ImageProducer)
