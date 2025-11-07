#ifndef ANALYSIS_IMAGE_CXX
#define ANALYSIS_IMAGE_CXX

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>

#include "AnalysisToolBase.h"
#include "Products/ImageProducts.h"
#include "Common/PandoraUtilities.h"
#include "Common/ProxyTypes.h"
#include "Imaging/SemanticClassifier.h"

#include <TDirectoryFile.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace analysis {

class ImageAnalysis : public AnalysisToolBase {
public:
  explicit ImageAnalysis(fhicl::ParameterSet const &p);
  virtual ~ImageAnalysis() = default;

  void configure(const fhicl::ParameterSet &p) override;
  void analyseEvent(const art::Event &event, bool is_data) override {}
  void analyseSlice(const art::Event &event,
                    std::vector<common::ProxyPfpElem_t> &pfp_pxy_vec,
                    bool is_data, bool is_selected) override;
  void setBranches(TTree *_tree) override;
  void resetTTree(TTree *_tree) override;

private:
  bool fDebug = false;
  art::InputTag fPFPproducer;
  art::InputTag fCLSproducer;
  art::InputTag fSLCproducer;
  art::InputTag fHITproducer;
  art::InputTag fWIREproducer;
  art::InputTag fMCPproducer;
  art::InputTag fBKTproducer;
  art::InputTag fImagesSliceTag;
  art::InputTag fImagesEventTag;
  float _reco_neutrino_vertex_x;
  float _reco_neutrino_vertex_y;
  float _reco_neutrino_vertex_z;
  std::vector<float> _detector_image_u;
  std::vector<float> _detector_image_v;
  std::vector<float> _detector_image_w;
  std::vector<int> _semantic_image_u;
  std::vector<int> _semantic_image_v;
  std::vector<int> _semantic_image_w;
  std::vector<float> _event_detector_image_u;
  std::vector<float> _event_detector_image_v;
  std::vector<float> _event_detector_image_w;
  std::vector<int> _event_semantic_image_u;
  std::vector<int> _event_semantic_image_v;
  std::vector<int> _event_semantic_image_w;
  float _event_adc_u;
  float _event_adc_v;
  float _event_adc_w;
  std::vector<int> _slice_semantic_counts_u;
  std::vector<int> _slice_semantic_counts_v;
  std::vector<int> _slice_semantic_counts_w;
  std::vector<int> _event_semantic_counts_u;
  std::vector<int> _event_semantic_counts_v;
  std::vector<int> _event_semantic_counts_w;
  bool _is_vtx_in_image_u;
  bool _is_vtx_in_image_v;
  bool _is_vtx_in_image_w;

  static std::vector<int> countLabels(const std::vector<int> &labels, size_t nlabels);
  void printSummary(const art::Event &event, bool is_data) const;
};

ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet &pset) {
  this->configure(pset);
}

void ImageAnalysis::configure(const fhicl::ParameterSet &p) {
  fDebug = p.get<bool>("Debug", false);
  fPFPproducer = p.get<art::InputTag>("PFPproducer");
  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fSLCproducer = p.get<art::InputTag>("SLCproducer");
  fHITproducer = p.get<art::InputTag>("HITproducer");
  fWIREproducer = p.get<art::InputTag>("WIREproducer");
  fMCPproducer = p.get<art::InputTag>("MCPproducer");
  fBKTproducer = p.get<art::InputTag>("BKTproducer");

  fImagesSliceTag = p.get<art::InputTag>(
      "ImagesSliceTag", art::InputTag{"imageprod", "primary_slice"});
  fImagesEventTag = p.get<art::InputTag>(
      "ImagesEventTag", art::InputTag{"imageprod", "all_hits"});
}

void ImageAnalysis::setBranches(TTree *_tree) {
  _tree->Branch("reco_neutrino_vertex_x", &_reco_neutrino_vertex_x,
                "reco_neutrino_vertex_x/F");
  _tree->Branch("reco_neutrino_vertex_y", &_reco_neutrino_vertex_y,
                "reco_neutrino_vertex_y/F");
  _tree->Branch("reco_neutrino_vertex_z", &_reco_neutrino_vertex_z,
                "reco_neutrino_vertex_z/F");
  _tree->Branch("detector_image_u", &_detector_image_u);
  _tree->Branch("detector_image_v", &_detector_image_v);
  _tree->Branch("detector_image_w", &_detector_image_w);
  _tree->Branch("semantic_image_u", &_semantic_image_u);
  _tree->Branch("semantic_image_v", &_semantic_image_v);
  _tree->Branch("semantic_image_w", &_semantic_image_w);
  _tree->Branch("event_detector_image_u", &_event_detector_image_u);
  _tree->Branch("event_detector_image_v", &_event_detector_image_v);
  _tree->Branch("event_detector_image_w", &_event_detector_image_w);
  _tree->Branch("event_semantic_image_u", &_event_semantic_image_u);
  _tree->Branch("event_semantic_image_v", &_event_semantic_image_v);
  _tree->Branch("event_semantic_image_w", &_event_semantic_image_w);
  _tree->Branch("event_adc_u", &_event_adc_u, "event_adc_u/F");
  _tree->Branch("event_adc_v", &_event_adc_v, "event_adc_v/F");
  _tree->Branch("event_adc_w", &_event_adc_w, "event_adc_w/F");
  _tree->Branch("slice_semantic_counts_u", &_slice_semantic_counts_u);
  _tree->Branch("slice_semantic_counts_v", &_slice_semantic_counts_v);
  _tree->Branch("slice_semantic_counts_w", &_slice_semantic_counts_w);
  _tree->Branch("event_semantic_counts_u", &_event_semantic_counts_u);
  _tree->Branch("event_semantic_counts_v", &_event_semantic_counts_v);
  _tree->Branch("event_semantic_counts_w", &_event_semantic_counts_w);
  _tree->Branch("is_vtx_in_image_u", &_is_vtx_in_image_u,
                "is_vtx_in_image_u/O");
  _tree->Branch("is_vtx_in_image_v", &_is_vtx_in_image_v,
                "is_vtx_in_image_v/O");
  _tree->Branch("is_vtx_in_image_w", &_is_vtx_in_image_w,
                "is_vtx_in_image_w/O");
}

void ImageAnalysis::resetTTree(TTree *_tree) {
  _reco_neutrino_vertex_x = std::numeric_limits<float>::quiet_NaN();
  _reco_neutrino_vertex_y = std::numeric_limits<float>::quiet_NaN();
  _reco_neutrino_vertex_z = std::numeric_limits<float>::quiet_NaN();
  _detector_image_u.clear();
  _detector_image_v.clear();
  _detector_image_w.clear();
  _semantic_image_u.clear();
  _semantic_image_v.clear();
  _semantic_image_w.clear();
  _event_detector_image_u.clear();
  _event_detector_image_v.clear();
  _event_detector_image_w.clear();
  _event_semantic_image_u.clear();
  _event_semantic_image_v.clear();
  _event_semantic_image_w.clear();
  _event_adc_u = std::numeric_limits<float>::quiet_NaN();
  _event_adc_v = std::numeric_limits<float>::quiet_NaN();
  _event_adc_w = std::numeric_limits<float>::quiet_NaN();
  _slice_semantic_counts_u.clear();
  _slice_semantic_counts_v.clear();
  _slice_semantic_counts_w.clear();
  _event_semantic_counts_u.clear();
  _event_semantic_counts_v.clear();
  _event_semantic_counts_w.clear();
  _is_vtx_in_image_u = false;
  _is_vtx_in_image_v = false;
  _is_vtx_in_image_w = false;
}

std::vector<int> ImageAnalysis::countLabels(const std::vector<int> &labels,
                                            size_t nlabels) {
  std::vector<int> counts(nlabels, 0);
  for (int v : labels) {
    if (v >= 0 && static_cast<size_t>(v) < nlabels) {
      ++counts[v];
    }
  }
  return counts;
}

void ImageAnalysis::analyseSlice(
    const art::Event &event, std::vector<common::ProxyPfpElem_t> &pfp_pxy_vec,
    bool is_data, bool is_selected) {
  for (const auto &pfp : pfp_pxy_vec) {
    if (pfp->IsPrimary()) {
      auto vtx = pfp.get<recob::Vertex>();
      if (vtx.size() > 0) {
        const auto &pos = vtx.at(0)->position();
        _reco_neutrino_vertex_x = pos.X();
        _reco_neutrino_vertex_y = pos.Y();
        _reco_neutrino_vertex_z = pos.Z();
      }
      break;
    }
  }

  auto sliceH = event.getValidHandle<std::vector<image::ImageProduct>>(fImagesSliceTag);
  auto eventH = event.getValidHandle<std::vector<image::ImageProduct>>(fImagesEventTag);

  auto assignPlane = [&](const image::ImageProduct &img, bool slice) {
    std::vector<float> *det_slice = nullptr;
    std::vector<int> *sem_slice = nullptr;
    std::vector<float> *det_event = nullptr;
    std::vector<int> *sem_event = nullptr;
    if (img.view == static_cast<int>(geo::kU)) {
      det_slice = &_detector_image_u;
      sem_slice = &_semantic_image_u;
      det_event = &_event_detector_image_u;
      sem_event = &_event_semantic_image_u;
    } else if (img.view == static_cast<int>(geo::kV)) {
      det_slice = &_detector_image_v;
      sem_slice = &_semantic_image_v;
      det_event = &_event_detector_image_v;
      sem_event = &_event_semantic_image_v;
    } else if (img.view == static_cast<int>(geo::kW)) {
      det_slice = &_detector_image_w;
      sem_slice = &_semantic_image_w;
      det_event = &_event_detector_image_w;
      sem_event = &_event_semantic_image_w;
    } else {
      return;
    }
    if (slice) {
      if (det_slice) *det_slice = img.adc;
      if (sem_slice) sem_slice->assign(img.semantic.begin(), img.semantic.end());
    } else {
      if (det_event) *det_event = img.adc;
      if (sem_event) sem_event->assign(img.semantic.begin(), img.semantic.end());
    }
  };

  for (const auto &pi : *sliceH) assignPlane(pi, true);
  for (const auto &pi : *eventH) assignPlane(pi, false);

  _event_adc_u = std::accumulate(_event_detector_image_u.begin(),
                                 _event_detector_image_u.end(), 0.0f);
  _event_adc_v = std::accumulate(_event_detector_image_v.begin(),
                                 _event_detector_image_v.end(), 0.0f);
  _event_adc_w = std::accumulate(_event_detector_image_w.begin(),
                                 _event_detector_image_w.end(), 0.0f);

  if (!is_data) {
    size_t nlabels = image::sem::SemanticClassifier::semantic_label_names.size();
    _slice_semantic_counts_u = countLabels(_semantic_image_u, nlabels);
    _slice_semantic_counts_v = countLabels(_semantic_image_v, nlabels);
    _slice_semantic_counts_w = countLabels(_semantic_image_w, nlabels);
    _event_semantic_counts_u = countLabels(_event_semantic_image_u, nlabels);
    _event_semantic_counts_v = countLabels(_event_semantic_image_v, nlabels);
    _event_semantic_counts_w = countLabels(_event_semantic_image_w, nlabels);
  }

  if (!std::isnan(_reco_neutrino_vertex_x)) {
    TVector3 vtx_pos_3d(_reco_neutrino_vertex_x, _reco_neutrino_vertex_y,
                        _reco_neutrino_vertex_z);
    TVector3 vtx_proj_u = common::ProjectToWireView(
        vtx_pos_3d.X(), vtx_pos_3d.Y(), vtx_pos_3d.Z(), common::TPC_VIEW_U);
    TVector3 vtx_proj_v = common::ProjectToWireView(
        vtx_pos_3d.X(), vtx_pos_3d.Y(), vtx_pos_3d.Z(), common::TPC_VIEW_V);
    TVector3 vtx_proj_w = common::ProjectToWireView(
        vtx_pos_3d.X(), vtx_pos_3d.Y(), vtx_pos_3d.Z(), common::TPC_VIEW_W);
    auto in_img = [](const image::ImageProduct &im, double drift, double wire) {
      bool in_row = (drift >= im.origin_y) &&
                    (drift < im.origin_y + im.pixel_h * static_cast<double>(im.height));
      bool in_col = (wire >= im.origin_x) &&
                    (wire < im.origin_x + im.pixel_w * static_cast<double>(im.width));
      return in_row && in_col;
    };
    const image::ImageProduct *U = nullptr;
    const image::ImageProduct *V = nullptr;
    const image::ImageProduct *W = nullptr;
    for (const auto &im : *sliceH) {
      if (im.view == static_cast<int>(geo::kU)) U = &im;
      if (im.view == static_cast<int>(geo::kV)) V = &im;
      if (im.view == static_cast<int>(geo::kW)) W = &im;
    }
    _is_vtx_in_image_u = (U && in_img(*U, vtx_proj_u.X(), vtx_proj_u.Z()));
    _is_vtx_in_image_v = (V && in_img(*V, vtx_proj_v.X(), vtx_proj_v.Z()));
    _is_vtx_in_image_w = (W && in_img(*W, vtx_proj_w.X(), vtx_proj_w.Z()));
  }

  if (fDebug) {
    printSummary(event, is_data);
  }
}

void ImageAnalysis::printSummary(const art::Event &event, bool is_data) const {
  const char *cat = "ImageAnalysis";
  const auto eid = event.id();
  mf::LogInfo(cat) << "=== ImageAnalysis summary: Run " << eid.run()
                   << " SubRun " << eid.subRun()
                   << " Event " << eid.event() << " ===";

  auto checkPair = [&](const char *name, const std::vector<float> &adc,
                       const std::vector<int> &sem) {
    bool mismatch = (!adc.empty() && !sem.empty() && adc.size() != sem.size());
    mf::LogInfo(cat) << "  " << name << " ADC=" << adc.size()
                     << "px, SEM=" << sem.size()
                     << "px" << (mismatch ? "  [SIZE-MISMATCH]" : "");
  };
  auto yesno = [](bool v) { return v ? "Y" : "N"; };

  checkPair("Slice U", _detector_image_u, _semantic_image_u);
  checkPair("Slice V", _detector_image_v, _semantic_image_v);
  checkPair("Slice W", _detector_image_w, _semantic_image_w);
  checkPair("Event U", _event_detector_image_u, _event_semantic_image_u);
  checkPair("Event V", _event_detector_image_v, _event_semantic_image_v);
  checkPair("Event W", _event_detector_image_w, _event_semantic_image_w);

  mf::LogInfo(cat) << "  Event ADC sums [U,V,W] = [" << _event_adc_u << ", "
                   << _event_adc_v << ", " << _event_adc_w << "]";

  mf::LogInfo(cat) << "  Reco vtx: (" << _reco_neutrino_vertex_x << ", "
                   << _reco_neutrino_vertex_y << ", "
                   << _reco_neutrino_vertex_z << ")";
  mf::LogInfo(cat) << "  Vertex in image [U,V,W] = ["
                   << yesno(_is_vtx_in_image_u) << ", "
                   << yesno(_is_vtx_in_image_v) << ", "
                   << yesno(_is_vtx_in_image_w) << "]";

  if (!is_data) {
    const auto &names = image::sem::SemanticClassifier::semantic_label_names;
    mf::LogInfo(cat) << "  Slice semantic counts by label:";
    for (size_t i = 0; i < names.size(); ++i) {
      int u = (i < _slice_semantic_counts_u.size()) ? _slice_semantic_counts_u[i] : 0;
      int v = (i < _slice_semantic_counts_v.size()) ? _slice_semantic_counts_v[i] : 0;
      int w = (i < _slice_semantic_counts_w.size()) ? _slice_semantic_counts_w[i] : 0;
      mf::LogInfo(cat) << "    [" << i << "] " << names[i]
                       << ": U=" << u << " V=" << v << " W=" << w;
    }
    mf::LogInfo(cat) << "  Event semantic counts by label:";
    for (size_t i = 0; i < names.size(); ++i) {
      int u = (i < _event_semantic_counts_u.size()) ? _event_semantic_counts_u[i] : 0;
      int v = (i < _event_semantic_counts_v.size()) ? _event_semantic_counts_v[i] : 0;
      int w = (i < _event_semantic_counts_w.size()) ? _event_semantic_counts_w[i] : 0;
      mf::LogInfo(cat) << "    [" << i << "] " << names[i]
                       << ": U=" << u << " V=" << v << " W=" << w;
    }
  }

  mf::LogInfo(cat) << "=== end ImageAnalysis summary ===";
}

DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif
