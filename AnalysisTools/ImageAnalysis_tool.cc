#ifndef ANALYSIS_IMAGE_CXX
#define ANALYSIS_IMAGE_CXX

#include "AnalysisToolBase.h"
#include "Common/Image.h"
#include "Common/ImageAlgorithm.h"
#include "Common/NpyUtils.h"
#include "Common/Pandora.h"
#include "Common/Types.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
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
#include "lardata/Utilities/FindManyInChainP.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <algorithm>
#include <array>
#include <cetlib/search_path.h>
#include <cmath>
#include <cstdlib>
#include <experimental/filesystem>
#include <fstream>
#include <iomanip>
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>
#include <limits.h>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::experimental::filesystem;

namespace analysis {

static inline bool is_remote_url(const std::string &s) {
  std::string_view sv{s};
  return sv.starts_with("http:") || sv.starts_with("https:");
}

static inline bool is_pnfs(const std::string &s) {
  std::string_view sv{s};
  return sv.starts_with("/pnfs/") || sv.starts_with("/eos/") ||
         sv.starts_with("/stash/");
}

    static std::optional<std::string>
    find_file_nearby(const fs::path &start, const std::string &filename,
                     int max_depth = 3) {
      fs::recursive_directory_iterator it(
          start, fs::directory_options::skip_permission_denied),
          end;
      for (; it != end; ++it) {
        if (it.depth() > max_depth) {
          it.pop();
          continue;
        }

        if (!fs::is_regular_file(it->path()))
          continue;
        if (it->path().filename() == filename)
          return it->path().string();
      }
      return std::nullopt;
    }

    static std::string resolve_weights(const std::string &weights_file,
                                       const std::string &work_dir,
                                       const std::string &base_dir,
                                       const std::string &scratch_dir,
                                       const std::string &model_name,
                                       bool fetch_ifdh) {
      fs::path wf(weights_file);

      if (wf.is_absolute() && fs::exists(wf))
        return wf.string();

      cet::search_path sp("FW_SEARCH_PATH");
      std::string resolved;
      if (sp.find_file(wf.string(), resolved))
        return resolved;

      if (!base_dir.empty()) {
        fs::path p = fs::path(base_dir) / wf;
        if (fs::exists(p))
          return p.string();
      }

      for (fs::path p : {fs::path(work_dir) / "weights" / wf.filename(),
                         fs::path(work_dir) / wf,
                         fs::path(scratch_dir) / "weights" / wf.filename()}) {
        if (fs::exists(p))
          return p.string();
      }

      if (auto hit = find_file_nearby(work_dir, wf.filename().string(), 3))
        return *hit;

      if (fetch_ifdh && is_pnfs(weights_file)) {
        fs::path dst =
            fs::path(scratch_dir) / (std::string("weights_") + model_name +
                                     "_" + wf.filename().string());
        std::ostringstream cmd;
        cmd << "ifdh cp -D " << weights_file << " " << dst.string();
        int rc = std::system(cmd.str().c_str());
        if (rc == 0 && fs::exists(dst))
          return dst.string();
      }

      std::ostringstream msg;
      msg << "Weights not found for model '" << model_name << "'\n"
          << "  requested: " << weights_file << "\n"
          << "  tried base: " << (fs::path(base_dir) / wf) << "\n"
          << "  tried cwd : " << (fs::path(work_dir) / wf)
          << " and weights/ subdir\n"
          << "  scratch   : "
          << (fs::path(scratch_dir) / "weights" / wf.filename()) << "\n"
          << "  (set WeightsBaseDir or WEIGHTS_BASE_DIR, or provide absolute "
             "path, or enable FetchWeightsWithIFDH and use PNFS path)";
      throw art::Exception(art::errors::Configuration) << msg.str();
    }

    struct ModelConfig {
      std::string name;
      std::string weights_file;
      std::string arch;
    };

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
      art::InputTag fSLCprodcuer;
      art::InputTag fHITproducer;
      art::InputTag fWIREproducer;
      art::InputTag fMCPproducer;
      art::InputTag fBKTproducer;
      std::string fBadChannelFile;
      std::set<unsigned int> fBadChannels;
      std::vector<ModelConfig> fModels;
      std::string fWeightsBaseDir;
      std::vector<std::string> fActiveModels;
      bool fFetchWeightsWithIFDH;
      int _image_width;
      int _image_height;
      float _adc_image_threshold;
      const geo::GeometryCore *_geo;
      const detinfo::DetectorProperties *_detp;
      float _drift_step;
      float _wire_pitch_u;
      float _wire_pitch_v;
      float _wire_pitch_w;
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
      std::unique_ptr<TruthLabelClassifier> _semantic_classifier;
      std::unordered_map<std::string, float> _inference_scores;

      void loadBadChannels(const std::string &filename);
      std::vector<art::Ptr<recob::Hit>> collectAllHits(const art::Event &event);
      std::vector<art::Ptr<recob::Hit>>
      collectSliceHits(const art::Event &event,
                       const std::vector<common::ProxyPfpElem_t> &pfp_pxy_vec);
      std::pair<double, double>
      calculateChargeCentroid(const art::Event &event, common::PandoraView view,
                              const std::vector<art::Ptr<recob::Hit>> &hits);
      void fillDetectorImage(
          const recob::Wire &wire, size_t wire_idx,
          const std::vector<ImageProperties> &properties,
          const art::FindManyP<recob::Hit> &wire_hit_assoc,
          const std::set<art::Ptr<recob::Hit>> &hit_set,
          const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
              &mcp_bkth_assoc,
          const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
          const std::map<int, size_t> &trackid_to_index,
          const std::vector<TruthLabelClassifier::TruthPrimaryLabel>
              &semantic_label_vector,
          std::vector<Image<float>> &detector_images,
          std::vector<Image<int>> &semantic_images, bool is_data,
          bool has_mcps);
      TruthLabelClassifier::TruthPrimaryLabel labelSemanticPixels(
          const art::Ptr<recob::Hit> &matched_hit,
          const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
              &mcp_bkth_assoc,
          const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
          const std::map<int, size_t> &trackid_to_index,
          const std::vector<TruthLabelClassifier::TruthPrimaryLabel>
              &semantic_label_vector,
          bool has_mcps) const;
      void constructPixelImages(const art::Event &event,
                                const std::vector<art::Ptr<recob::Hit>> &hits,
                                const std::vector<ImageProperties> &properties,
                                std::vector<Image<float>> &detector_images,
                                std::vector<Image<int>> &semantic_images,
                                bool is_data);
      float runInference(const std::vector<Image<float>> &detector_images,
                         const std::string &absolute_scratch_dir,
                         const std::string &work_dir, const std::string &arch,
                         const std::string &weights_file);
    };

    ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet &pset) {
      this->configure(pset);
    }

    void ImageAnalysis::configure(const fhicl::ParameterSet &p) {
      fPFPproducer = p.get<art::InputTag>("PFPproducer");
      fCLSproducer = p.get<art::InputTag>("CLSproducer");
      fSLCprodcuer = p.get<art::InputTag>("SLCproducer");
      fHITproducer = p.get<art::InputTag>("HITproducer");
      fWIREproducer = p.get<art::InputTag>("WIREproducer");
      fMCPproducer = p.get<art::InputTag>("MCPproducer");
      fBKTproducer = p.get<art::InputTag>("BKTproducer");
      fBadChannelFile = p.get<std::string>("BadChannelFile");
      if (!fBadChannelFile.empty()) {
        cet::search_path sp("FW_SEARCH_PATH");
        std::string resolved_path;
        if (sp.find_file(fBadChannelFile, resolved_path)) {
          fBadChannelFile = resolved_path;
          mf::LogInfo("ImageAnalysis")
              << "Bad channel file resolved to '" << fBadChannelFile << "'";
        } else {
          mf::LogWarning("ImageAnalysis")
              << "Bad channel file '" << fBadChannelFile
              << "' not found via FW_SEARCH_PATH";
          fBadChannelFile.clear();
        }
        char cwd[PATH_MAX];
        if (getcwd(cwd, sizeof(cwd))) {
          mf::LogInfo log("ImageAnalysis");
          log << "Job CWD: " << cwd << '\n';
          try {
            for (const auto &entry : fs::directory_iterator(cwd)) {
              log << "  " << entry.path().filename().string() << '\n';
            }
          } catch (const fs::filesystem_error &e) {
            mf::LogWarning("ImageAnalysis")
                << "Unable to list directory contents: " << e.what();
          }
        }
      }
      fWeightsBaseDir = p.get<std::string>("WeightsBaseDir", "");
      fActiveModels = p.get<std::vector<std::string>>("ActiveModels", {});
      fFetchWeightsWithIFDH = p.get<bool>("FetchWeightsWithIFDH", false);
      if (const char *wbd = std::getenv("WEIGHTS_BASE_DIR"))
        fWeightsBaseDir = wbd;
      if (const char *am = std::getenv("ACTIVE_MODELS")) {
        fActiveModels.clear();
        std::stringstream ss(am);
        std::string tok;
        while (std::getline(ss, tok, ','))
          if (!tok.empty())
            fActiveModels.push_back(tok);
      }
      auto model_psets = p.get<std::vector<fhicl::ParameterSet>>("Models", {});
      for (const auto &ps : model_psets) {
        fModels.push_back({ps.get<std::string>("name"),
                           ps.get<std::string>("weights_file"),
                           ps.get<std::string>("arch")});
      }
      _image_width = 512;
      _image_height = 512;
      _adc_image_threshold = 1.0;
      _geo = art::ServiceHandle<geo::Geometry>()->provider();
      _detp =
          art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
      auto clock =
          art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
      double tick_period = clock->TPCClock().TickPeriod();
      double drift_velocity = _detp->DriftVelocity();
      _drift_step = tick_period * drift_velocity * 1e1;
      _wire_pitch_u = _geo->WirePitch(geo::kU);
      _wire_pitch_v = _geo->WirePitch(geo::kV);
      _wire_pitch_w = _geo->WirePitch(geo::kW);
      _semantic_classifier =
          std::make_unique<TruthLabelClassifier>(fMCPproducer);
      if (!fBadChannelFile.empty()) {
        this->loadBadChannels(fBadChannelFile);
      }
    }

    void ImageAnalysis::loadBadChannels(const std::string &filename) {
      fBadChannels.clear();
      std::ifstream infile(filename);
      if (!infile.is_open()) {
        throw art::Exception(art::errors::Configuration)
            << "Could not open bad channels file: " << filename;
      }
      std::string line;
      while (std::getline(infile, line)) {
        if (line.empty() || line.find('#') != std::string::npos) {
          continue;
        }
        std::stringstream ss(line);
        unsigned int first_ch, second_ch;
        ss >> first_ch;
        if (ss >> second_ch) {
          for (unsigned int ch = first_ch; ch <= second_ch; ++ch) {
            fBadChannels.insert(ch);
          }
        } else {
          fBadChannels.insert(first_ch);
        }
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
      for (const auto &model : fModels) {
        _tree->Branch(("inference_score_" + model.name).c_str(),
                      &_inference_scores[model.name], "score/F");
      }
    }

    void ImageAnalysis::resetTTree(TTree *_tree) {
      _reco_neutrino_vertex_x = std::numeric_limits<float>::lowest();
      _reco_neutrino_vertex_y = std::numeric_limits<float>::lowest();
      _reco_neutrino_vertex_z = std::numeric_limits<float>::lowest();
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
      _event_adc_u = std::numeric_limits<float>::lowest();
      _event_adc_v = std::numeric_limits<float>::lowest();
      _event_adc_w = std::numeric_limits<float>::lowest();
      _slice_semantic_counts_u.clear();
      _slice_semantic_counts_v.clear();
      _slice_semantic_counts_w.clear();
      _event_semantic_counts_u.clear();
      _event_semantic_counts_v.clear();
      _event_semantic_counts_w.clear();
      _is_vtx_in_image_u = false;
      _is_vtx_in_image_v = false;
      _is_vtx_in_image_w = false;
      for (auto &kv : _inference_scores) {
        kv.second = std::numeric_limits<float>::quiet_NaN();
      }
    }

    void ImageAnalysis::analyseSlice(
        const art::Event &event,
        std::vector<common::ProxyPfpElem_t> &pfp_pxy_vec, bool is_data,
        bool is_selected) {
      for (const auto &pfp : pfp_pxy_vec) {
        if (pfp->IsPrimary()) {
          const auto &vtx = pfp.get<recob::Vertex>();
          if (vtx.size() > 0) {
            const auto &vtx_pos = vtx[0]->position();
            _reco_neutrino_vertex_x = vtx_pos.X();
            _reco_neutrino_vertex_y = vtx_pos.Y();
            _reco_neutrino_vertex_z = vtx_pos.Z();
          }
          break;
        }
      }

      std::vector<art::Ptr<recob::Hit>> all_hits = this->collectAllHits(event);
      std::vector<art::Ptr<recob::Hit>> neutrino_hits =
          this->collectSliceHits(event, pfp_pxy_vec);

      auto [centroid_wire_u, centroid_drift_u] = this->calculateChargeCentroid(
          event, common::TPC_VIEW_U, neutrino_hits);
      auto [centroid_wire_v, centroid_drift_v] = this->calculateChargeCentroid(
          event, common::TPC_VIEW_V, neutrino_hits);
      auto [centroid_wire_w, centroid_drift_w] = this->calculateChargeCentroid(
          event, common::TPC_VIEW_W, neutrino_hits);

      std::vector<ImageProperties> properties;
      properties.emplace_back(centroid_wire_u, centroid_drift_u, _image_width,
                              _image_height, _drift_step, _wire_pitch_u,
                              geo::kU);
      properties.emplace_back(centroid_wire_v, centroid_drift_v, _image_width,
                              _image_height, _drift_step, _wire_pitch_v,
                              geo::kV);
      properties.emplace_back(centroid_wire_w, centroid_drift_w, _image_width,
                              _image_height, _drift_step, _wire_pitch_w,
                              geo::kW);

      std::vector<Image<float>> detector_images;
      std::vector<Image<int>> semantic_images;
      this->constructPixelImages(event, neutrino_hits, properties,
                                 detector_images, semantic_images, is_data);

      std::vector<Image<float>> event_detector_images;
      std::vector<Image<int>> event_semantic_images;
      this->constructPixelImages(event, all_hits, properties,
                                 event_detector_images, event_semantic_images,
                                 is_data);

      _detector_image_u = detector_images[0].data();
      _detector_image_v = detector_images[1].data();
      _detector_image_w = detector_images[2].data();
      _semantic_image_u = semantic_images[0].data();
      _semantic_image_v = semantic_images[1].data();
      _semantic_image_w = semantic_images[2].data();
      _event_detector_image_u = event_detector_images[0].data();
      _event_detector_image_v = event_detector_images[1].data();
      _event_detector_image_w = event_detector_images[2].data();
      _event_semantic_image_u = event_semantic_images[0].data();
      _event_semantic_image_v = event_semantic_images[1].data();
      _event_semantic_image_w = event_semantic_images[2].data();
      _event_adc_u = std::accumulate(_event_detector_image_u.begin(),
                                     _event_detector_image_u.end(), 0.0f);
      _event_adc_v = std::accumulate(_event_detector_image_v.begin(),
                                     _event_detector_image_v.end(), 0.0f);
      _event_adc_w = std::accumulate(_event_detector_image_w.begin(),
                                     _event_detector_image_w.end(), 0.0f);

      if (!is_data) {
        for (size_t i = 0;
             i < TruthLabelClassifier::truth_primary_label_names.size(); ++i) {
          _slice_semantic_counts_u.push_back(std::count(
              _semantic_image_u.begin(), _semantic_image_u.end(), i));
          _slice_semantic_counts_v.push_back(std::count(
              _semantic_image_v.begin(), _semantic_image_v.end(), i));
          _slice_semantic_counts_w.push_back(std::count(
              _semantic_image_w.begin(), _semantic_image_w.end(), i));
          _event_semantic_counts_u.push_back(
              std::count(_event_semantic_image_u.begin(),
                         _event_semantic_image_u.end(), i));
          _event_semantic_counts_v.push_back(
              std::count(_event_semantic_image_v.begin(),
                         _event_semantic_image_v.end(), i));
          _event_semantic_counts_w.push_back(
              std::count(_event_semantic_image_w.begin(),
                         _event_semantic_image_w.end(), i));
        }
      }

      std::string scratch_dir = ".";
      const char *scratch_dir_cstr = std::getenv("_CONDOR_SCRATCH_DIR");
      if (scratch_dir_cstr) {
        scratch_dir = std::string(scratch_dir_cstr);
      }
      char scratch_path_buffer[PATH_MAX];
      realpath(scratch_dir.c_str(), scratch_path_buffer);
      std::string absolute_scratch_dir(scratch_path_buffer);

      char cwd_buffer[PATH_MAX];
      getcwd(cwd_buffer, sizeof(cwd_buffer));
      std::string work_dir(cwd_buffer);

      std::vector<ModelConfig> todo;
      if (fActiveModels.empty())
        todo = fModels;
      else {
        std::set<std::string> want(fActiveModels.begin(), fActiveModels.end());
        for (auto const &m : fModels)
          if (want.count(m.name))
            todo.push_back(m);
      }

      for (auto const &m : todo) {
        std::string wpath = resolve_weights(
            m.weights_file, work_dir, fWeightsBaseDir, absolute_scratch_dir,
            m.name, fFetchWeightsWithIFDH);
        mf::LogInfo("ImageAnalysis")
            << "Model " << m.name << " using weights: " << wpath;
        float score = this->runInference(detector_images, absolute_scratch_dir,
                                         work_dir, m.arch, wpath);
        _inference_scores[m.name] = score;
      }

      if (_reco_neutrino_vertex_x != std::numeric_limits<float>::lowest()) {
        TVector3 vtx_pos_3d(_reco_neutrino_vertex_x, _reco_neutrino_vertex_y,
                            _reco_neutrino_vertex_z);
        TVector3 vtx_proj_u = common::ProjectToWireView(
            vtx_pos_3d.X(), vtx_pos_3d.Y(), vtx_pos_3d.Z(), common::TPC_VIEW_U);
        _is_vtx_in_image_u =
            (properties[0].row(vtx_proj_u.X()) != static_cast<size_t>(-1) &&
             properties[0].col(vtx_proj_u.Z()) != static_cast<size_t>(-1));
        TVector3 vtx_proj_v = common::ProjectToWireView(
            vtx_pos_3d.X(), vtx_pos_3d.Y(), vtx_pos_3d.Z(), common::TPC_VIEW_V);
        _is_vtx_in_image_v =
            (properties[1].row(vtx_proj_v.X()) != static_cast<size_t>(-1) &&
             properties[1].col(vtx_proj_v.Z()) != static_cast<size_t>(-1));
        TVector3 vtx_proj_w = common::ProjectToWireView(
            vtx_pos_3d.X(), vtx_pos_3d.Y(), vtx_pos_3d.Z(), common::TPC_VIEW_W);
        _is_vtx_in_image_w =
            (properties[2].row(vtx_proj_w.X()) != static_cast<size_t>(-1) &&
             properties[2].col(vtx_proj_w.Z()) != static_cast<size_t>(-1));
      }
    }

    std::vector<art::Ptr<recob::Hit>>
    ImageAnalysis::collectAllHits(const art::Event &event) {
      std::vector<art::Ptr<recob::Hit>> all_hits;
      auto hit_handle =
          event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
      for (size_t i = 0; i < hit_handle->size(); ++i) {
        all_hits.emplace_back(hit_handle, i);
      }
      return all_hits;
    }

    std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectSliceHits(
        const art::Event &event,
        const std::vector<common::ProxyPfpElem_t> &pfp_pxy_vec) {
      std::vector<art::Ptr<recob::Hit>> neutrino_hits;
      auto pfpHandle =
          event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
      art::FindManyP<recob::Slice> sliceAssoc(pfpHandle, event, fPFPproducer);
      size_t pfpIndex = pfp_pxy_vec[0].index();
      const auto &slices = sliceAssoc.at(pfpIndex);
      if (slices.empty())
        return {};
      const art::Ptr<recob::Slice> &slice = slices[0];
      auto sliceHandle =
          event.getValidHandle<std::vector<recob::Slice>>(fSLCprodcuer);
      art::FindManyP<recob::Hit> hitAssoc(sliceHandle, event, fSLCprodcuer);
      const std::vector<art::Ptr<recob::Hit>> &sliceHits =
          hitAssoc.at(slice.key());
      neutrino_hits.reserve(sliceHits.size());
      for (const auto &hit : sliceHits) {
        neutrino_hits.push_back(hit);
      }
      return neutrino_hits;
    }

    std::pair<double, double> ImageAnalysis::calculateChargeCentroid(
        const art::Event &event, common::PandoraView view,
        const std::vector<art::Ptr<recob::Hit>> &hits) {
      double sum_charge = 0.0;
      double sum_wire = 0.0;
      double sum_drift = 0.0;
      for (const auto &hit : hits) {
        if (fBadChannels.count(hit->Channel())) {
          continue;
        }
        if (common::GetPandoraView(hit) != view)
          continue;
        double charge = hit->Integral();
        TVector3 hit_pos = common::GetPandoraHitPosition(event, hit, view);
        sum_charge += charge;
        sum_wire += hit_pos.Z() * charge;
        sum_drift += hit_pos.X() * charge;
      }
      if (sum_charge == 0.0)
        return {0.0, 0.0};
      return {sum_wire / sum_charge, sum_drift / sum_charge};
    }

    TruthLabelClassifier::TruthPrimaryLabel
    ImageAnalysis::labelSemanticPixels(
        const art::Ptr<recob::Hit> &matched_hit,
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
            &mcp_bkth_assoc,
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
        const std::map<int, size_t> &trackid_to_index,
        const std::vector<TruthLabelClassifier::TruthPrimaryLabel>
            &semantic_label_vector,
        bool has_mcps) const {
      TruthLabelClassifier::TruthPrimaryLabel semantic_pixel_label =
          TruthLabelClassifier::TruthPrimaryLabel::Cosmic;
      if (!has_mcps || !mcp_vector.isValid())
        return semantic_pixel_label;
      std::vector<art::Ptr<simb::MCParticle>> mcp_particles_ass_to_hit;
      std::vector<anab::BackTrackerHitMatchingData const *> bkth_data_ass_to_hit;
      mcp_bkth_assoc.get(matched_hit.key(), mcp_particles_ass_to_hit,
                         bkth_data_ass_to_hit);
      if (!bkth_data_ass_to_hit.empty()) {
        float max_ide_fraction = -1.0;
        int best_match_track_id = -1;
        for (size_t i_bkth = 0; i_bkth < bkth_data_ass_to_hit.size(); ++i_bkth) {
          if (bkth_data_ass_to_hit[i_bkth] &&
              bkth_data_ass_to_hit[i_bkth]->ideFraction > max_ide_fraction) {
            max_ide_fraction = bkth_data_ass_to_hit[i_bkth]->ideFraction;
            best_match_track_id = mcp_particles_ass_to_hit[i_bkth]->TrackId();
          }
        }
        if (best_match_track_id != -1) {
          auto it_trackid = trackid_to_index.find(best_match_track_id);
          if (it_trackid != trackid_to_index.end()) {
            size_t particle_mcp_idx = it_trackid->second;
            if (particle_mcp_idx < semantic_label_vector.size()) {
              semantic_pixel_label =
                  semantic_label_vector[particle_mcp_idx];
            }
          }
        }
      }
      return semantic_pixel_label;
    }

    void ImageAnalysis::fillDetectorImage(
        const recob::Wire &wire, size_t wire_idx,
        const std::vector<ImageProperties> &properties,
        const art::FindManyP<recob::Hit> &wire_hit_assoc,
        const std::set<art::Ptr<recob::Hit>> &hit_set,
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
            &mcp_bkth_assoc,
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
        const std::map<int, size_t> &trackid_to_index,
        const std::vector<TruthLabelClassifier::TruthPrimaryLabel>
            &semantic_label_vector,
        std::vector<Image<float>> &detector_images,
        std::vector<Image<int>> &semantic_images, bool is_data,
        bool has_mcps) {
      auto ch_id = wire.Channel();
      if (fBadChannels.count(ch_id)) {
        return;
      }
      std::vector<geo::WireID> wire_ids = _geo->ChannelToWire(ch_id);
      if (wire_ids.empty())
        return;
      geo::View_t view = _geo->View(wire_ids.front().planeID());
      size_t view_idx = static_cast<size_t>(view);
      const geo::WireGeo *wire_geo = _geo->WirePtr(wire_ids.front());
      TVector3 center = wire_geo->GetCenter();
      TVector3 wire_center(center.X(), center.Y(), center.Z());
      double wire_coord =
          (view == geo::kW)
              ? wire_center.Z()
              : (view == geo::kU)
                    ? (wire_center.Z() * std::cos(1.04719758034) -
                       wire_center.Y() * std::sin(1.04719758034))
                    : (wire_center.Z() * std::cos(-1.04719758034) -
                       wire_center.Y() * std::sin(-1.04719758034));
      auto hits_for_wire = wire_hit_assoc.at(wire_idx);
      std::vector<art::Ptr<recob::Hit>> filtered_hits;
      for (const auto &hit : hits_for_wire) {
        if (hit_set.count(hit)) {
          filtered_hits.push_back(hit);
        }
      }
      std::sort(filtered_hits.begin(), filtered_hits.end(),
                [](const art::Ptr<recob::Hit> &a,
                   const art::Ptr<recob::Hit> &b) {
                  return a->StartTick() < b->StartTick();
                });
      size_t hit_index = 0;
      for (const auto &range : wire.SignalROI().get_ranges()) {
        const auto &adcs = range.data();
        int start_tick = range.begin_index();
        for (size_t adc_index = 0; adc_index < adcs.size(); ++adc_index) {
          int tick = start_tick + adc_index;
          double x_drift_pos = _detp->ConvertTicksToX(
              static_cast<double>(tick), wire_ids.front().planeID());
          size_t row = properties[view_idx].row(x_drift_pos);
          size_t col = properties[view_idx].col(wire_coord);
          if (row == static_cast<size_t>(-1) ||
              col == static_cast<size_t>(-1))
            continue;
          while (hit_index < filtered_hits.size() &&
                 filtered_hits[hit_index]->EndTick() <= tick) {
            ++hit_index;
          }
          if (hit_index < filtered_hits.size() &&
              filtered_hits[hit_index]->StartTick() <= tick &&
              tick < filtered_hits[hit_index]->EndTick()) {
            if (adcs[adc_index] > _adc_image_threshold) {
              detector_images[view_idx].set(row, col, adcs[adc_index]);
              if (!is_data) {
                TruthLabelClassifier::TruthPrimaryLabel semantic_pixel_label =
                    labelSemanticPixels(filtered_hits[hit_index], mcp_bkth_assoc,
                                        mcp_vector, trackid_to_index,
                                        semantic_label_vector, has_mcps);
                semantic_images[view_idx].set(
                    row, col, static_cast<int>(semantic_pixel_label), false);
              }
            }
          }
        }
      }
    }

    void ImageAnalysis::constructPixelImages(
        const art::Event &event, const std::vector<art::Ptr<recob::Hit>> &hits,
        const std::vector<ImageProperties> &properties,
        std::vector<Image<float>> &detector_images,
        std::vector<Image<int>> &semantic_images, bool is_data) {
      detector_images.clear();
      semantic_images.clear();
      for (const auto &prop : properties) {
        Image<float> detector_image(prop);
        detector_image.clear(0.0);
        detector_images.push_back(std::move(detector_image));
        Image<int> semantic_image(prop);
        semantic_image.clear(
            static_cast<int>(TruthLabelClassifier::TruthPrimaryLabel::Empty));
        semantic_images.push_back(std::move(semantic_image));
      }

      auto wire_vector =
          event.getValidHandle<std::vector<recob::Wire>>(fWIREproducer);
      auto hit_vector =
          event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
      art::Handle<std::vector<simb::MCParticle>> mcp_vector;
      bool has_mcps = event.getByLabel(fMCPproducer, mcp_vector);
      art::FindManyP<recob::Hit> wire_hit_assoc(wire_vector, event,
                                                fHITproducer);
      art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
          mcp_bkth_assoc(hit_vector, event, fBKTproducer);

      std::vector<TruthLabelClassifier::TruthPrimaryLabel>
          semantic_label_vector;
      if (!is_data && has_mcps && mcp_vector.isValid() &&
          _semantic_classifier) {
        semantic_label_vector = _semantic_classifier->classifyParticles(event);
      }

      std::map<int, size_t> trackid_to_index;
      if (!is_data && has_mcps && mcp_vector.isValid()) {
        for (size_t i = 0; i < mcp_vector->size(); ++i) {
          trackid_to_index[mcp_vector->at(i).TrackId()] = i;
        }
      }

      std::set<art::Ptr<recob::Hit>> hit_set(hits.begin(), hits.end());
      for (size_t wire_idx = 0; wire_idx < wire_vector->size(); ++wire_idx) {
        fillDetectorImage(wire_vector->at(wire_idx), wire_idx, properties,
                          wire_hit_assoc, hit_set, mcp_bkth_assoc, mcp_vector,
                          trackid_to_index, semantic_label_vector,
                          detector_images, semantic_images, is_data, has_mcps);
      }
    }

    float ImageAnalysis::runInference(
        const std::vector<Image<float>> &detector_images,
        const std::string &absolute_scratch_dir, const std::string &work_dir,
        const std::string &arch, const std::string &weights_file) {
      using std::string;
      string npy_in = absolute_scratch_dir + "/detector_images.npy";
      string temp_out = absolute_scratch_dir + "/temp_test_out.txt";
      string script_stdout = absolute_scratch_dir + "/py_script.out";
      string script_stderr = absolute_scratch_dir + "/py_script.err";

      std::vector<std::vector<float>> images;
      images.push_back(detector_images[0].data());
      images.push_back(detector_images[1].data());
      images.push_back(detector_images[2].data());
      common::save_npy_f32_2d(npy_in, images);

      string container = "/cvmfs/uboone.opensciencegrid.org/containers/"
                         "lantern_v2_me_06_03_prod";
      cet::search_path sp("FW_SEARCH_PATH");
      string wrapper_script;
      if (!sp.find_file("run_strangeness_inference.sh", wrapper_script)) {
        mf::LogWarning("ImageAnalysis")
            << "run_strangeness_inference.sh not found via FW_SEARCH_PATH";
        return 0.f;
      }

      std::string bind_paths = absolute_scratch_dir;
      bool need_cvmfs = true;
      if (const char *bind_env = std::getenv("APPTAINER_BINDPATH")) {
        if (std::string(bind_env).find("/cvmfs") != std::string::npos)
          need_cvmfs = false;
      } else {
        std::ifstream mounts("/proc/mounts");
        std::string line;
        while (std::getline(mounts, line)) {
          if (line.find(" /cvmfs ") != std::string::npos) {
            need_cvmfs = false;
            break;
          }
        }
      }
      if (need_cvmfs)
        bind_paths = "/cvmfs," + bind_paths;

      std::ostringstream cmd;
      cmd << "apptainer exec --cleanenv --bind " << bind_paths << " "
          << container << " "
          << "/bin/bash " << wrapper_script << " "
          << "--npy " << npy_in << " "
          << "--output " << temp_out << " "
          << "--arch " << arch << " "
          << "--weights " << weights_file << " > " << script_stdout << " 2> "
          << script_stderr;

      mf::LogInfo("ImageAnalysis") << "Executing inference: " << cmd.str();

      auto start = std::chrono::steady_clock::now();
      int code = std::system(cmd.str().c_str());
      auto end = std::chrono::steady_clock::now();
      double duration = std::chrono::duration<double>(end - start).count();

      bool success = (code == 0) && fs::exists(temp_out);
      if (!success) {
        std::ifstream error_stream(script_stderr);
        std::string error_message(
            (std::istreambuf_iterator<char>(error_stream)), {});
        throw art::Exception(art::errors::LogicError)
            << "Inference script failed with exit code " << code
            << "\n--- Script stderr ---\n"
            << error_message << "\n--- End Script stderr ---";
      }

      std::ifstream result_stream(temp_out.c_str());
      if (!result_stream) {
        throw art::Exception(art::errors::LogicError)
            << "Could not open temporary result file: " << temp_out;
      }
      float score;
      result_stream >> score;

      mf::LogInfo("ImageAnalysis")
          << "Inference time: " << duration << " seconds";
      mf::LogInfo("ImageAnalysis") << "Predicted score: " << score;

      std::error_code ec;
      fs::remove(npy_in, ec);
      fs::remove(temp_out, ec);
      fs::remove(script_stdout, ec);
      fs::remove(script_stderr, ec);

      return score;
    }

    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
    } // namespace analysis

#endif
