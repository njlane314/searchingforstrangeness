#ifndef ANALYSIS_IMAGE_CXX
#define ANALYSIS_IMAGE_CXX

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
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <cetlib/search_path.h>
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>

#include "AnalysisToolBase.h"
#include "Imaging/Image.h"
#include "Imaging/SemanticPixelClassifier.h"
#include "Imaging/ImageAlgo.h"
#include "Common/NpyUtils.h"
#include "Common/PandoraUtilities.h"
#include "Common/ProxyTypes.h"

#include <TDirectoryFile.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits.h>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

namespace analysis {

namespace {
inline bool fileExists(const std::string &p) {
    struct stat sb;
    return ::stat(p.c_str(), &sb) == 0 && S_ISREG(sb.st_mode);
}
inline bool dirExists(const std::string &p) {
    struct stat sb;
    return ::stat(p.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode);
}
inline bool isAbs(const std::string &p) {
    return !p.empty() && p.front() == '/';
}

static double geometryInformedChargeFraction(int width_px, int height_px) {
    const double Dx = std::hypot(width_px / 2.0, height_px / 2.0);
    const double R_px = 0.5 * Dx;
    const double lambda = std::min(1.0, R_px / Dx);
    return std::min(1.0, 0.5 * M_PI * lambda * lambda);
}

static std::pair<double, double> centroidClosestRCharge(
    const art::Event &event, common::PandoraView view,
    const std::vector<art::Ptr<recob::Hit>> &hits, double r_keep,
    const std::set<unsigned int> &bad_channels) {
    struct P {
        double Z, X, q;
    };
    std::vector<P> pts;
    pts.reserve(hits.size());
    std::vector<std::pair<double, double>> Zw, Xw;
    Zw.reserve(hits.size());
    Xw.reserve(hits.size());

    for (auto const &h : hits) {
        if (bad_channels.count(h->Channel()))
            continue;
        if (common::GetPandoraView(h) != view)
            continue;
        double q = std::max(0.f, h->Integral());
        if (q <= 0)
            continue;
        TVector3 p = common::GetPandoraHitPosition(event, h, view);
        pts.push_back({p.Z(), p.X(), q});
        Zw.emplace_back(p.Z(), q);
        Xw.emplace_back(p.X(), q);
    }
    if (pts.empty())
        return {0.0, 0.0};

    auto wmedian = [](std::vector<std::pair<double, double>> v) {
        std::sort(v.begin(), v.end(),
                  [](auto &a, auto &b) { return a.first < b.first; });
        double W = 0;
        for (auto &t : v)
            W += t.second;
        double acc = 0;
        for (auto &t : v) {
            acc += t.second;
            if (acc >= 0.5 * W)
                return t.first;
        }
        return v.back().first;
    };

    const double z0 = wmedian(Zw), x0 = wmedian(Xw);

    struct DW {
        double d, w;
    };
    std::vector<DW> dw;
    dw.reserve(pts.size());
    double Wtot = 0;
    for (auto const &p : pts) {
        double d = std::hypot(p.Z - z0, p.X - x0);
        dw.push_back({d, p.q});
        Wtot += p.q;
    }

    std::sort(dw.begin(), dw.end(),
              [](auto &a, auto &b) { return a.d < b.d; });
    double target = r_keep * Wtot, acc = 0;
    size_t k = 0;
    for (; k < dw.size() && acc < target; ++k)
        acc += dw[k].w;

    double W = 0, Zs = 0, Xs = 0;
    const double d_cut = (k == 0 ? 0.0 : dw[k - 1].d);
    for (auto const &p : pts) {
        double d = std::hypot(p.Z - z0, p.X - x0);
        if (d <= d_cut) {
            W += p.q;
            Zs += p.q * p.Z;
            Xs += p.q * p.X;
        }
    }
    if (W == 0)
        return {z0, x0};
    return {Zs / W, Xs / W};
}

} // namespace

class ImageAnalysis : public AnalysisToolBase {
  public:
    explicit ImageAnalysis(fhicl::ParameterSet const &p);
    virtual ~ImageAnalysis() = default;

    void configure(const fhicl::ParameterSet &p) override;
    void analyseEvent(const art::Event &event, bool is_data) override {}
    void analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &pfp_pxy_vec, bool is_data, bool is_selected) override;
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
    std::string fBadChannelFile;
    std::set<unsigned int> fBadChannels;
    std::vector<ModelConfig> fModels;
    std::string fWeightsBaseDir;
    std::vector<std::string> fActiveModels;
    std::string fAssetsBaseDir;
    std::string fInferenceWrapper;
    std::string fWorkDir;
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
    std::unique_ptr<SemanticPixelClassifier> _semantic_classifier;
    std::unique_ptr<ImageAlgo> _image_algo;
    std::unordered_map<std::string, float> _inference_scores;

    void loadBadChannels(const std::string &filename);
    std::vector<art::Ptr<recob::Hit>> collectAllHits(const art::Event &event);
    std::vector<art::Ptr<recob::Hit>> collectSliceHits(const art::Event &event, const std::vector<common::ProxyPfpElem_t> &pfp_pxy_vec);
    std::string resolveAssetPath(const std::string &relOrAbs, bool must_exist = true) const;
};

ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet &pset) {
    this->configure(pset);
}

void ImageAnalysis::configure(const fhicl::ParameterSet &p) {
    char cwd_buffer[PATH_MAX];
    getcwd(cwd_buffer, sizeof(cwd_buffer));
    fWorkDir = std::string(cwd_buffer);

    fAssetsBaseDir = p.get<std::string>("AssetsBaseDir", "");
    if (const char *env = std::getenv("ASSETS_BASE_DIR"))
        fAssetsBaseDir = env;
    if (!fAssetsBaseDir.empty() && !isAbs(fAssetsBaseDir))
        fAssetsBaseDir = joinPath(fWorkDir, fAssetsBaseDir);

    std::string badRel = p.get<std::string>("BadChannelFile", "badchannels.txt");
    if (const char *env = std::getenv("IA_BADCHANNELS"))
        badRel = env;
    fBadChannelFile = resolveAssetPath(badRel, true);
    mf::LogInfo("ImageAnalysis") << "Bad channel file: " << fBadChannelFile;
    this->loadBadChannels(fBadChannelFile);

    fWeightsBaseDir = p.get<std::string>("WeightsBaseDir", "weights");
    if (const char *env = std::getenv("WEIGHTS_BASE_DIR"))
        fWeightsBaseDir = env;
    if (!isAbs(fWeightsBaseDir))
        fWeightsBaseDir = resolveAssetPath(fWeightsBaseDir, false);
    if (!dirExists(fWeightsBaseDir)) {
        throw art::Exception(art::errors::Configuration)
            << "Weights base dir '" << fWeightsBaseDir << "' does not exist.";
    }
    mf::LogInfo("ImageAnalysis") << "Weights base dir: " << fWeightsBaseDir;

    std::string wrapperRel = p.get<std::string>("InferenceWrapper", "scripts/inference_wrapper.sh");
    if (const char *env = std::getenv("IA_INFERENCE_WRAPPER"))
        wrapperRel = env;
    fInferenceWrapper = resolveAssetPath(wrapperRel, true);
    mf::LogInfo("ImageAnalysis") << "Inference wrapper: " << fInferenceWrapper;

    fPFPproducer = p.get<art::InputTag>("PFPproducer");
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fSLCproducer = p.get<art::InputTag>("SLCproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fWIREproducer = p.get<art::InputTag>("WIREproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fBKTproducer = p.get<art::InputTag>("BKTproducer");

    fActiveModels = p.get<std::vector<std::string>>("ActiveModels", {});
    if (const char *am = std::getenv("ACTIVE_MODELS")) {
        fActiveModels.clear();
        std::stringstream ss(am);
        std::string tok;
        while (std::getline(ss, tok, ','))
            if (!tok.empty())
                fActiveModels.push_back(tok);
    }

    auto model_psets = p.get<std::vector<fhicl::ParameterSet>>("Models", {});
    fModels.clear();
    for (const auto &ps : model_psets) {
        fModels.push_back({ps.get<std::string>("name"),
                           ps.get<std::string>("weights_file"),
                           ps.get<std::string>("arch")});
    }

    _image_width = 512;
    _image_height = 512;
    // Set the default ADC threshold for image production
    _adc_image_threshold = 4.0;
    _geo = art::ServiceHandle<geo::Geometry>()->provider();
    _detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
    double tick_period = clock->TPCClock().TickPeriod();
    double drift_velocity = _detp->DriftVelocity();
    _drift_step = tick_period * drift_velocity * 1e1;
    _wire_pitch_u = _geo->WirePitch(geo::kU);
    _wire_pitch_v = _geo->WirePitch(geo::kV);
    _wire_pitch_w = _geo->WirePitch(geo::kW);
    _semantic_classifier = std::make_unique<SemanticPixelClassifier>(fMCPproducer);
    _image_algo = std::make_unique<ImageAlgo>(fWIREproducer, fHITproducer,
                                             fMCPproducer, fBKTproducer,
                                             _adc_image_threshold, fWeightsBaseDir,
                                             fInferenceWrapper, fAssetsBaseDir,
                                             fModels, fActiveModels, _geo, _detp,
                                             fWorkDir);
}

std::string ImageAnalysis::resolveAssetPath(const std::string &relOrAbs, bool must_exist) const {
    if (relOrAbs.empty()) return {};
    if (isAbs(relOrAbs) && (!must_exist || fileExists(relOrAbs) || dirExists(relOrAbs)))
        return relOrAbs;
    if (!fAssetsBaseDir.empty()) {
        auto cand = joinPath(fAssetsBaseDir, relOrAbs);
        if (!must_exist || fileExists(cand) || dirExists(cand))
            return cand;
    }
    if (!fWorkDir.empty()) {
        auto cand = joinPath(fWorkDir, relOrAbs);
        if (!must_exist || fileExists(cand) || dirExists(cand))
            return cand;
    }
    cet::search_path sp("FW_SEARCH_PATH");
    std::string resolved;
    if (sp.find_file(relOrAbs, resolved))
        return resolved;
    if (must_exist) {
        throw art::Exception(art::errors::Configuration)
            << "Could not resolve path '" << relOrAbs
            << "'. Tried under AssetsBaseDir='" << fAssetsBaseDir
            << "', CWD='" << fWorkDir << "', and FW_SEARCH_PATH.";
    }
    return relOrAbs;
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
    for (auto &kv : _inference_scores) {
        kv.second = std::numeric_limits<float>::quiet_NaN();
    }
}

void ImageAnalysis::analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &pfp_pxy_vec, bool is_data, bool is_selected) {
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

    const double r_keep = geometryInformedChargeFraction(
        _image_width, _image_height);

    auto [centroid_wire_u, centroid_drift_u] = centroidClosestRCharge(
        event, common::TPC_VIEW_U, neutrino_hits, r_keep, fBadChannels);
    auto [centroid_wire_v, centroid_drift_v] = centroidClosestRCharge(
        event, common::TPC_VIEW_V, neutrino_hits, r_keep, fBadChannels);
    auto [centroid_wire_w, centroid_drift_w] = centroidClosestRCharge(
        event, common::TPC_VIEW_W, neutrino_hits, r_keep, fBadChannels);

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
    _image_algo->produceImages(event, neutrino_hits, properties, is_data,
                              _semantic_classifier.get(), fBadChannels,
                              detector_images, semantic_images);

    std::vector<Image<float>> event_detector_images;
    std::vector<Image<int>> event_semantic_images;
    _image_algo->produceImages(event, all_hits, properties, is_data,
                              _semantic_classifier.get(), fBadChannels,
                              event_detector_images, event_semantic_images);

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
             i < SemanticPixelClassifier::semantic_label_names.size(); ++i) {
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

    auto scores = _image_algo->runInference(detector_images, absolute_scratch_dir);
    for (auto const &kv : scores) {
        _inference_scores[kv.first] = kv.second;
    }

    if (!std::isnan(_reco_neutrino_vertex_x)) {
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

std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectAllHits(const art::Event &event) {
    std::vector<art::Ptr<recob::Hit>> all_hits;
    auto hit_handle =
        event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    for (size_t i = 0; i < hit_handle->size(); ++i) {
        all_hits.emplace_back(hit_handle, i);
    }
    return all_hits;
}

std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectSliceHits(const art::Event &event, const std::vector<common::ProxyPfpElem_t> &pfp_pxy_vec) {
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
        event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    art::FindManyP<recob::Hit> hitAssoc(sliceHandle, event, fSLCproducer);
    const std::vector<art::Ptr<recob::Hit>> &sliceHits =
        hitAssoc.at(slice.key());
    neutrino_hits.reserve(sliceHits.size());
    for (const auto &hit : sliceHits) {
        neutrino_hits.push_back(hit);
    }
    return neutrino_hits;
}





DEFINE_ART_CLASS_TOOL(ImageAnalysis)
} // namespace analysis

#endif
