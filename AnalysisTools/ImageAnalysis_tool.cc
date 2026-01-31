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
#include "Products/InferencePerf.h"
#include "Products/InferencePred.h"
#include "Products/ImageProducts.h"
#include "Common/PandoraUtilities.h"
#include "Common/ProxyTypes.h"
#include "Imaging/SemanticClassifier.h"

#include <TDirectoryFile.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <algorithm>
#include <cstdint>
#include <cmath>
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
  art::InputTag fPFPproducer;
  art::InputTag fCLSproducer;
  art::InputTag fSLCproducer;
  art::InputTag fHITproducer;
  art::InputTag fWIREproducer;
  art::InputTag fMCPproducer;
  art::InputTag fBKTproducer;
  art::InputTag fImagesSliceTag;
  art::InputTag fInferencePredTag;
  art::InputTag fInferencePerfTag;
  float _reco_neutrino_vertex_x;
  float _reco_neutrino_vertex_y;
  float _reco_neutrino_vertex_z;
  std::vector<float> _detector_image_u;
  std::vector<float> _detector_image_v;
  std::vector<float> _detector_image_w;
  std::vector<int32_t> _semantic_image_u;
  std::vector<int32_t> _semantic_image_v;
  std::vector<int32_t> _semantic_image_w;
  std::vector<int> _slice_semantic_counts_u;
  std::vector<int> _slice_semantic_counts_v;
  std::vector<int> _slice_semantic_counts_w;
  bool _is_vtx_in_image_u;
  bool _is_vtx_in_image_v;
  bool _is_vtx_in_image_w;

  // Inference output: predictions
  std::vector<std::string> _inf_model;
  std::vector<int> _inf_n_scores;
  std::vector<float> _inf_scores;

  // Inference output: performance metrics
  std::vector<float> _inf_t_write_req_ms;
  std::vector<float> _inf_t_exec_total_ms;
  std::vector<float> _inf_t_read_resp_ms;
  std::vector<float> _inf_t_child_total_ms;
  std::vector<float> _inf_t_child_setup_ms;
  std::vector<float> _inf_t_child_infer_ms;
  std::vector<float> _inf_t_child_post_ms;
  std::vector<float> _inf_child_max_rss_mb;

  static std::vector<int> countLabels(const std::vector<int32_t> &labels,
                                      size_t nlabels);
};

ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet &pset) {
  this->configure(pset);
}

void ImageAnalysis::configure(const fhicl::ParameterSet &p) {
  fPFPproducer = p.get<art::InputTag>("PFPproducer");
  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fSLCproducer = p.get<art::InputTag>("SLCproducer");
  fHITproducer = p.get<art::InputTag>("HITproducer");
  fWIREproducer = p.get<art::InputTag>("WIREproducer");
  fMCPproducer = p.get<art::InputTag>("MCPproducer");
  fBKTproducer = p.get<art::InputTag>("BKTproducer");

  fImagesSliceTag =
      p.get<art::InputTag>("ImagesSliceTag",
                           art::InputTag{"imageprod", "NuSlice"});

  if (p.has_key("InferencePredTag")) {
    fInferencePredTag = p.get<art::InputTag>("InferencePredTag");
  }

  if (p.has_key("InferencePerfTag")) {
    fInferencePerfTag = p.get<art::InputTag>("InferencePerfTag");
  }
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
  _tree->Branch("slice_semantic_counts_u", &_slice_semantic_counts_u);
  _tree->Branch("slice_semantic_counts_v", &_slice_semantic_counts_v);
  _tree->Branch("slice_semantic_counts_w", &_slice_semantic_counts_w);
  _tree->Branch("is_vtx_in_image_u", &_is_vtx_in_image_u,
                "is_vtx_in_image_u/O");
  _tree->Branch("is_vtx_in_image_v", &_is_vtx_in_image_v,
                "is_vtx_in_image_v/O");
  _tree->Branch("is_vtx_in_image_w", &_is_vtx_in_image_w,
                "is_vtx_in_image_w/O");

  _tree->Branch("inf_model", &_inf_model);
  _tree->Branch("inf_n_scores", &_inf_n_scores);
  _tree->Branch("inf_scores", &_inf_scores);

  _tree->Branch("inf_t_write_req_ms", &_inf_t_write_req_ms);
  _tree->Branch("inf_t_exec_total_ms", &_inf_t_exec_total_ms);
  _tree->Branch("inf_t_read_resp_ms", &_inf_t_read_resp_ms);
  _tree->Branch("inf_t_child_total_ms", &_inf_t_child_total_ms);
  _tree->Branch("inf_t_child_setup_ms", &_inf_t_child_setup_ms);
  _tree->Branch("inf_t_child_infer_ms", &_inf_t_child_infer_ms);
  _tree->Branch("inf_t_child_post_ms", &_inf_t_child_post_ms);
  _tree->Branch("inf_child_max_rss_mb", &_inf_child_max_rss_mb);
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
  _slice_semantic_counts_u.clear();
  _slice_semantic_counts_v.clear();
  _slice_semantic_counts_w.clear();
  _is_vtx_in_image_u = false;
  _is_vtx_in_image_v = false;
  _is_vtx_in_image_w = false;

  _inf_model.clear();
  _inf_n_scores.clear();
  _inf_scores.clear();

  _inf_t_write_req_ms.clear();
  _inf_t_exec_total_ms.clear();
  _inf_t_read_resp_ms.clear();
  _inf_t_child_total_ms.clear();
  _inf_t_child_setup_ms.clear();
  _inf_t_child_infer_ms.clear();
  _inf_t_child_post_ms.clear();
  _inf_child_max_rss_mb.clear();
}

std::vector<int> ImageAnalysis::countLabels(const std::vector<int32_t> &labels,
                                            size_t nlabels) {
  std::vector<int> counts(nlabels, 0);
  for (int32_t v : labels) {
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

  auto assignPlane = [&](const image::ImageProduct &img, bool slice) {
    std::vector<float> *det_slice = nullptr;
    std::vector<int32_t> *sem_slice = nullptr;
    if (img.view == static_cast<int>(geo::kU)) {
      det_slice = &_detector_image_u;
      sem_slice = &_semantic_image_u;
    } else if (img.view == static_cast<int>(geo::kV)) {
      det_slice = &_detector_image_v;
      sem_slice = &_semantic_image_v;
    } else if (img.view == static_cast<int>(geo::kW)) {
      det_slice = &_detector_image_w;
      sem_slice = &_semantic_image_w;
    } else {
      return;
    }
    if (slice) {
      if (det_slice) *det_slice = img.adc;
      if (sem_slice) sem_slice->assign(img.semantic.begin(), img.semantic.end());
    }
  };

  for (const auto &pi : *sliceH) assignPlane(pi, true);

  if (!is_data) {
    size_t nlabels = image::SemanticClassifier::semantic_label_names.size();
    _slice_semantic_counts_u = countLabels(_semantic_image_u, nlabels);
    _slice_semantic_counts_v = countLabels(_semantic_image_v, nlabels);
    _slice_semantic_counts_w = countLabels(_semantic_image_w, nlabels);
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

  if (!fInferencePredTag.label().empty()) {
    art::Handle<image::InferencePredProduct> predH;
    event.getByLabel(fInferencePredTag, predH);
    if (predH.isValid()) {
      for (auto const &mp : predH->per_model) {
        _inf_model.push_back(mp.model);
        _inf_n_scores.push_back(static_cast<int>(mp.scores.size()));
        _inf_scores.insert(_inf_scores.end(), mp.scores.begin(), mp.scores.end());
      }
    } else {
      mf::LogWarning("ImageAnalysis")
          << "InferencePredProduct with tag \"" << fInferencePredTag.encode()
          << "\" not found or invalid.";
    }
  }

  if (!fInferencePerfTag.label().empty()) {
    art::Handle<image::InferencePerfProduct> perfH;
    event.getByLabel(fInferencePerfTag, perfH);
    if (perfH.isValid()) {
      for (auto const &mp : perfH->per_model) {
        _inf_t_write_req_ms.push_back(mp.t_write_req_ms);
        _inf_t_exec_total_ms.push_back(mp.t_exec_total_ms);
        _inf_t_read_resp_ms.push_back(mp.t_read_resp_ms);
        _inf_t_child_total_ms.push_back(mp.t_child_total_ms);
        _inf_t_child_setup_ms.push_back(mp.t_child_setup_ms);
        _inf_t_child_infer_ms.push_back(mp.t_child_infer_ms);
        _inf_t_child_post_ms.push_back(mp.t_child_post_ms);
        _inf_child_max_rss_mb.push_back(mp.child_max_rss_mb);
      }
    } else {
      mf::LogWarning("ImageAnalysis")
          << "InferencePerfProduct with tag \"" << fInferencePerfTag.encode()
          << "\" not found or invalid.";
    }
  }

}

DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif
