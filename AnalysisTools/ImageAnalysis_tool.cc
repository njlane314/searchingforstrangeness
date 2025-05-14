#include <vector>
#include <string>
#include <map>
#include <array>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <limits>

#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TVector3.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ zowelplugin/tool_macros.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>
#include "lardataobj/AnalysisBase/MVAOutput.h"


#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/Utilities/FindManyInChainP.h"


#include "../CommonDefs/Pandora.h"
#include "../CommonDefs/Types.h"

namespace analysis {

    constexpr double SIXTY_DEGREES_RAD = 1.04719758034;

    class ImageProperties {
    public:
        ImageProperties() = default;
        ImageProperties(double center_x_cm, double center_y_cm, size_t width_pixels, size_t height_pixels, double pixel_height_cm, double pixel_width_cm, geo::View_t view)
            : center_x_cm_(center_x_cm), center_y_cm_(center_y_cm),
              height_pixels_(height_pixels), width_pixels_(width_pixels),
              pixel_width_cm_(pixel_width_cm), pixel_height_cm_(pixel_height_cm), view_(view) {
            origin_x_cm_ = center_x_cm_ - (static_cast<double>(width_pixels_) * pixel_width_cm_) / 2.0;
            origin_y_cm_ = center_y_cm_ - (static_cast<double>(height_pixels_) * pixel_height_cm_) / 2.0;
        }

        size_t index(size_t row, size_t col) const {
            if (row >= height_pixels_ || col >= width_pixels_) return static_cast<size_t>(-1);
            return col * height_pixels_ + row;
        }

        size_t col(double x_cm) const {
            if (x_cm < origin_x_cm_ || x_cm >= origin_x_cm_ + static_cast<double>(width_pixels_) * pixel_width_cm_) return static_cast<size_t>(-1);
            return static_cast<size_t>((x_cm - origin_x_cm_) / pixel_width_cm_);
        }

        size_t row(double y_cm) const {
            if (y_cm < origin_y_cm_ || y_cm >= origin_y_cm_ + static_cast<double>(height_pixels_) * pixel_height_cm_) return static_cast<size_t>(-1);
            return static_cast<size_t>((y_cm - origin_y_cm_) / pixel_height_cm_);
        }

        size_t height() const { return height_pixels_; }
        size_t width() const { return width_pixels_; }
        double pixelWidthCm() const { return pixel_width_cm_; }
        double pixelHeightCm() const { return pixel_height_cm_; }
        geo::View_t view() const { return view_; }
        double originXCm() const { return origin_x_cm_; }
        double originYCm() const { return origin_y_cm_; }
        double maxXCm() const { return origin_x_cm_ + static_cast<double>(width_pixels_) * pixel_width_cm_; }
        double maxYCm() const { return origin_y_cm_ + static_cast<double>(height_pixels_) * pixel_height_cm_; }

    private:
        double center_x_cm_, center_y_cm_;
        double origin_x_cm_, origin_y_cm_;
        size_t height_pixels_, width_pixels_;
        double pixel_width_cm_, pixel_height_cm_;
        geo::View_t view_ {geo::kUnknown};
    };

    template <typename T>
    class Image {
    public:
        Image() = default;
        Image(const ImageProperties& prop)
            : prop_(prop), pixels_(prop.height() * prop.width(), T(0)) {}

        void set(size_t row, size_t col, T value, bool accumulate = true) {
            size_t idx = prop_.index(row, col);
            if (idx != static_cast<size_t>(-1)) {
                if (accumulate)
                    pixels_[idx] += value;
                else
                    pixels_[idx] = value;
            }
        }

        T get(size_t row, size_t col) const {
            size_t idx = prop_.index(row, col);
            return (idx != static_cast<size_t>(-1)) ? pixels_[idx] : T(0);
        }

        void clear(T value = T(0)) {
            std::fill(pixels_.begin(), pixels_.end(), value);
        }

        const std::vector<T>& data() const {
            return pixels_;
        }

        geo::View_t view() const { return prop_.view(); }
        size_t height() const { return prop_.height(); }
        size_t width() const { return prop_.width(); }
        const ImageProperties& properties() const { return prop_; }

    private:
        ImageProperties prop_;
        std::vector<T> pixels_;
    };

    enum class TruthPrimaryLabel {
        Empty = 0, Cosmic, Electron, Muon,
        ChargedPion, NeutralPion, Proton,
        ChargedKaon, NeutralKaon, Lambda, ChargedSigma,
        OtherHadron, OtherParticle
    };

    const std::array<std::string, 13> truth_primary_label_names = {
        "Empty", "Cosmic", "Electron", "Muon",
        "ChargedPion", "NeutralPion", "Proton",
        "ChargedKaon", "NeutralKaon", "Lambda", "ChargedSigma",
        "OtherHadron", "OtherParticle"
    };

    inline TruthPrimaryLabel getTruthPrimaryLabelFromPDG(int pdg_code_signed) {
        int pdg_code = std::abs(pdg_code_signed);

        if (pdg_code == 11) return TruthPrimaryLabel::Electron;
        if (pdg_code == 13) return TruthPrimaryLabel::Muon;
        if (std::abs(pdg_code_signed) == 211) return TruthPrimaryLabel::ChargedPion;
        if (pdg_code == 111) return TruthPrimaryLabel::NeutralPion;
        if (pdg_code == 2212) return TruthPrimaryLabel::Proton;
        if (pdg_code == 321) return TruthPrimaryLabel::ChargedKaon
        if (pdg_code == 311) return TruthPrimaryLabel::NeutralKaon;
        if (pdg_code == 3122) return TruthPrimaryLabel::Lambda;
        if (pdg_code == 3222 || pdg_code == 3112) return TruthPrimaryLabel::ChargedSigma;
        
        if (pdg_code > 1000 && pdg_code < 10000 && (pdg_code / 1000 % 2 != 0 || pdg_code / 100 % 2 != 0 || pdg_code / 10 % 2 != 0 ) ) {
             return TruthPrimaryLabel::OtherHadron;
        }
        return TruthPrimaryLabel::OtherParticle;
    }


    inline void processTruthPrimaryParticleRecursive(
        size_t particle_idx,
        const std::vector<simb::MCParticle>& mc_particles,
        std::vector<TruthPrimaryLabel>& particle_labels_vec,
        const std::unordered_map<int, size_t>& track_id_to_idx_map,
        TruthPrimaryLabel current_primary_label ) {
        if (particle_idx >= particle_labels_vec.size()) return;
        particle_labels_vec[particle_idx] = current_primary_label;
        const auto& mc_particle = mc_particles[particle_idx];
        for (int i = 0; i < mc_particle.NumberDaughters(); ++i) {
            if (auto it = track_id_to_idx_map.find(mc_particle.Daughter(i)); it != track_id_to_idx_map.end()) {
                if (it->second < mc_particles.size()) {
                     processTruthPrimaryParticleRecursive(it->second, mc_particles, particle_labels_vec, track_id_to_idx_map, current_primary_label);
                }
            }
        }
    }

    inline std::vector<TruthPrimaryLabel> classifyTruthPrimaryLabels(
        const art::Event& event, const art::InputTag& mc_particle_label_tag ) {
        auto mc_particle_handle = event.getValidHandle<std::vector<simb::MCParticle>>(mc_particle_label_tag);
        const auto& mc_particles = *mc_particle_handle;
        std::unordered_map<int, size_t> track_id_to_idx_map;
        for (size_t i = 0; i < mc_particles.size(); ++i) {
            track_id_to_idx_map[mc_particles[i].TrackId()] = i;
        }
        std::vector<TruthPrimaryLabel> particle_labels_vec(mc_particles.size(), TruthPrimaryLabel::Empty);
        for (size_t i = 0; i < mc_particles.size(); ++i) {
            if (mc_particles[i].Mother() == 0) {
                if (auto it = track_id_to_idx_map.find(mc_particles[i].TrackId()); it != track_id_to_idx_map.end()) {
                    TruthPrimaryLabel label = getTruthPrimaryLabelFromPDG(mc_particles[i].PdgCode());
                    if (mc_particles[i].Process() == "primaryCosmic" || (mc_particles[i].StatusCode() == 0 && mc_particles[i].Mother()==0 && mc_particles[i].Process().find("decay") != std::string::npos) ){
                        label = TruthPrimaryLabel::Cosmic;
                    }
                    processTruthPrimaryParticleRecursive(it->second, mc_particles, particle_labels_vec, track_id_to_idx_map, label);
                }
            }
        }
        return particle_labels_vec;
    }

    enum class RecoLabel {
        empty, cosmic, MIP, HIP, shower, michel, diffuse, invisible
    };

    const std::array<std::string, 8> reco_label_names = {
        "empty", "cosmic", "MIP", "HIP", "shower", "michel", "diffuse", "invisible"
    };

    inline std::pair<RecoLabel, RecoLabel> computeRecoLabel(
        const std::vector<simb::MCParticle>& mc_particles, const simb::MCParticle& current_mc_part,
        RecoLabel sl_from_parent, const std::map<int, size_t>& track_id_to_idx_map,
        double gamma_threshold_gev, double hadron_threshold_gev, double lepton_threshold_gev ) {

        if (sl_from_parent != RecoLabel::empty) {
            return {sl_from_parent, sl_from_parent};
        }

        RecoLabel sl = RecoLabel::invisible;
        RecoLabel slc = RecoLabel::empty;

        int pdg_code = current_mc_part.PdgCode();
        double momentum_gev = current_mc_part.P();
        const std::string& start_process = current_mc_part.Process();
        const std::string& end_process = current_mc_part.EndProcess();

        int parent_pdg_code = 0;
        if (current_mc_part.Mother() != 0) {
            auto it = track_id_to_idx_map.find(current_mc_part.Mother());
            if (it != track_id_to_idx_map.end()) {
                parent_pdg_code = mc_particles[it->second].PdgCode();
            }
        }

        if (std::abs(pdg_code) == 13) {
             if (momentum_gev > lepton_threshold_gev) sl = RecoLabel::MIP; else sl = RecoLabel::diffuse;
        } else if (std::abs(pdg_code) == 211) {
            sl = RecoLabel::MIP;
        } else if (std::abs(pdg_code) == 321 || (std::abs(pdg_code) == 2212 && momentum_gev >= hadron_threshold_gev)) {
            sl = RecoLabel::HIP;
        } else if (std::abs(pdg_code) == 11) {
            if (start_process == "primary" && momentum_gev >= lepton_threshold_gev) {
                sl = RecoLabel::shower;
                slc = RecoLabel::shower;
            } else if (std::abs(parent_pdg_code) == 13 && (start_process == "muMinusCaptureAtRest" || start_process == "muPlusCaptureAtRest" || start_process == "Decay")) {
                sl = RecoLabel::michel;
                slc = RecoLabel::michel;
            } else if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt") {
                if (momentum_gev >= gamma_threshold_gev) {
                    sl = RecoLabel::shower;
                    slc = RecoLabel::shower;
                } else {
                    sl = RecoLabel::diffuse;
                }
            } else {
                 if (momentum_gev >= lepton_threshold_gev) sl = RecoLabel::shower; else sl = RecoLabel::diffuse;
            }
        } else if (pdg_code == 22) {
            if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt" || start_process == "primary") {
                if (momentum_gev >= gamma_threshold_gev) {
                    sl = RecoLabel::shower;
                    slc = RecoLabel::shower;
                } else {
                    sl = RecoLabel::diffuse;
                }
            } else {
                sl = RecoLabel::diffuse;
            }
        } else if (std::abs(pdg_code) == 2212 && momentum_gev < hadron_threshold_gev) {
            sl = RecoLabel::diffuse;
        }
        return {sl, slc};
    }

    inline void processRecoLabelRecursive(
        size_t particle_idx, const std::vector<simb::MCParticle>& mc_particles,
        RecoLabel sl_from_parent, std::vector<RecoLabel>& particle_labels_vec,
        const std::map<int, size_t>& track_id_to_idx_map,
        double gamma_thresh, double hadron_thresh, double lepton_thresh) {
        const auto& mc_part = mc_particles[particle_idx];
        auto [sl, slc] = computeRecoLabel(mc_particles, mc_part, sl_from_parent, track_id_to_idx_map, gamma_thresh, hadron_thresh, lepton_thresh);
        particle_labels_vec[particle_idx] = sl;
        for (int i = 0; i < mc_part.NumberDaughters(); ++i) {
            int daughter_id = mc_part.Daughter(i);
            auto it = track_id_to_idx_map.find(daughter_id);
            if (it != track_id_to_idx_map.end()) {
                processRecoLabelRecursive(it->second, mc_particles, slc, particle_labels_vec, track_id_to_idx_map, gamma_thresh, hadron_thresh, lepton_thresh);
            }
        }
    }

    inline std::vector<RecoLabel> classifyRecoLabels(
        const art::Event& event, const art::InputTag& mc_particle_label_tag,
        double gamma_threshold_gev, double hadron_threshold_gev, double lepton_threshold_gev ) {
        auto mc_particle_handle = event.getValidHandle<std::vector<simb::MCParticle>>(mc_particle_label_tag);
        const auto& mc_particles = *mc_particle_handle;
        std::map<int, size_t> track_id_to_idx_map;
        for (size_t i = 0; i < mc_particles.size(); ++i) {
            track_id_to_idx_map[mc_particles[i].TrackId()] = i;
        }
        std::vector<RecoLabel> particle_labels_vec(mc_particles.size(), RecoLabel::empty);
        for (size_t i = 0; i < mc_particles.size(); ++i) {
            if (mc_particles[i].Mother() == 0) {
                processRecoLabelRecursive(i, mc_particles, RecoLabel::empty, particle_labels_vec, track_id_to_idx_map, gamma_threshold_gev, hadron_threshold_gev, lepton_threshold_gev);
            }
        }
        return particle_labels_vec;
    }

    template <typename T_score, typename A_score>
    int arg_max_score(std::vector<T_score, A_score> const& vec) {
        if (vec.empty()) return -1; 
        return static_cast<int>(std::distance(vec.begin(), std::max_element(vec.begin(), vec.end())));
    }


    class AnalysisToolBase {
    public:
        virtual ~AnalysisToolBase() = default;
        virtual void configure(const fhicl::ParameterSet& p) = 0;
        virtual void analyseEvent(art::Event const& e, bool _is_data) = 0;
        virtual void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) = 0;
        virtual void setBranches(TTree* tree) = 0;
        virtual void resetTTree(TTree* tree) = 0;
    };

    class ImageAnalysis : public AnalysisToolBase {
    public:
        explicit ImageAnalysis(fhicl::ParameterSet const& p);
        void configure(const fhicl::ParameterSet& p) override;
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) override;
        void setBranches(TTree* tree) override;

        void analyseEvent(art::Event const& e, bool _is_data) override {}
        void resetTTree(TTree* tree) override;

    private:
        art::InputTag fPFPproducer;
        art::InputTag fCLSproducer;
        art::InputTag fHITproducer;
        art::InputTag fWIREproducer;
        art::InputTag fMCPproducer;
        std::string fBackTrackerLabel;
        double fGammaThresholdGeV;
        double fHadronThresholdGeV;
        double fLeptonThresholdGeV;
        bool fProcessMC;

        int _image_width_pixels;
        int _image_height_pixels;
        float _adc_image_threshold;

        const geo::GeometryCore* _geo_serv;
        const detinfo::DetectorProperties* _det_prop_serv;

        float _drift_step_cm_per_tick;
        float _wire_pitch_u_cm;
        float _wire_pitch_v_cm;
        float _wire_pitch_w_cm;

        TTree* _image_tree_ptr;

        int img_event_id;
        int img_slice_id;
        int img_nrows;
        int img_ncols;
        bool img_is_vertex_in_image;

        std::vector<std::vector<float>> img_raw_adc_images;
        std::vector<std::vector<int>> img_cheated_reco_label_images;
        std::vector<std::vector<int>> img_truth_primary_label_images;


        std::vector<art::Ptr<recob::Hit>> collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);
        std::pair<double, double> calculateChargeWeightedCentroid(const art::Event& e, common::PandoraView target_view, const std::vector<art::Ptr<recob::Hit>>& hits);
        void constructPixelImages(const art::Event& e, const std::vector<ImageProperties>& view_img_props,
                                 std::vector<Image<float>>& out_raw_adc_imgs,
                                 std::vector<Image<int>>& out_cheated_reco_imgs,
                                 std::vector<Image<int>>& out_truth_primary_imgs);
    };

    ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet& pset) :
        _geo_serv(nullptr), _det_prop_serv(nullptr), _image_tree_ptr(nullptr) {
        this->configure(pset);
    }

    void ImageAnalysis::configure(const fhicl::ParameterSet& p) {
        fPFPproducer = p.get<art::InputTag>("PFPproducer");
        fCLSproducer = p.get<art::InputTag>("CLSproducer");
        fHITproducer = p.get<art::InputTag>("HITproducer");
        fWIREproducer = p.get<art::InputTag>("WIREproducer");
        fMCPproducer = p.get<art::InputTag>("MCPproducer");
        fBackTrackerLabel = p.get<std::string>("BackTrackerLabel", "gaushit");
        fGammaThresholdGeV = p.get<double>("GammaThresholdGeV", 0.02);
        fHadronThresholdGeV = p.get<double>("HadronThresholdGeV", 0.1);
        fLeptonThresholdGeV = p.get<double>("LeptonThresholdGeV", 0.01);
        fProcessMC = p.get<bool>("ProcessMC", false);

        _image_width_pixels = p.get<int>("ImageWidthPixels", 512);
        _image_height_pixels = p.get<int>("ImageHeightPixels", 512);
        _adc_image_threshold = p.get<float>("ADCImageThreshold", 4.0);

        _geo_serv = art::ServiceHandle<geo::Geometry>()->provider();
        _det_prop_serv = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        auto const clock_serv = art::ServiceHandle<detinfo::DetectorClocksService const>()->provider();

        double tick_period_us = clock_serv->TPCClock().TickPeriod();
        double drift_velocity_cm_per_us = _det_prop_serv->DriftVelocity();
        _drift_step_cm_per_tick = tick_period_us * drift_velocity_cm_per_us;

        _wire_pitch_u_cm = _geo_serv->WirePitch(geo::kU);
        _wire_pitch_v_cm = _geo_serv->WirePitch(geo::kV);
        _wire_pitch_w_cm = _geo_serv->WirePitch(geo::kW);

        img_raw_adc_images.resize(3);
        if (fProcessMC) {
            img_cheated_reco_label_images.resize(3);
            img_truth_primary_label_images.resize(3);
        }
    }
    
    void ImageAnalysis::setBranches(TTree* tree) {
        _image_tree_ptr = tree;
        if (!_image_tree_ptr) return;
        _image_tree_ptr->Branch("img_event_id", &img_event_id, "img_event_id/I");
        _image_tree_ptr->Branch("img_slice_id", &img_slice_id, "img_slice_id/I");
        _image_tree_ptr->Branch("img_nrows", &img_nrows, "img_nrows/I");
        _image_tree_ptr->Branch("img_ncols", &img_ncols, "img_ncols/I");
        _image_tree_ptr->Branch("img_raw_adc_images", &img_raw_adc_images);
        _image_tree_ptr->Branch("img_is_vertex_in_image", &img_is_vertex_in_image, "img_is_vertex_in_image/O");

        if (fProcessMC) {
            _image_tree_ptr->Branch("img_cheated_reco_label_images", &img_cheated_reco_label_images);
            _image_tree_ptr->Branch("img_truth_primary_label_images", &img_truth_primary_label_images);
        }
    }
    
    void ImageAnalysis::resetTTree(TTree* tree) {
        img_event_id = std::numeric_limits<int>::lowest();
        img_slice_id = std::numeric_limits<int>::lowest();
        img_nrows = 0;
        img_ncols = 0;
        img_is_vertex_in_image = false;

        for(auto& v : img_raw_adc_images) v.clear();
        if (fProcessMC) {
            for(auto& v : img_cheated_reco_label_images) v.clear();
            for(auto& v : img_truth_primary_label_images) v.clear();
        }
    }


    void ImageAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) {
        img_event_id = e.id().event();
        img_slice_id = (slice_pfp_v.empty() || slice_pfp_v[0].isNull()) ? -1 : static_cast<int>(slice_pfp_v[0]->Parent());


        std::vector<art::Ptr<recob::Hit>> current_slice_hits = this->collectSliceHits(e, slice_pfp_v);

        auto [centroid_wire_u, centroid_drift_u] = this->calculateChargeWeightedCentroid(e, common::TPC_VIEW_U, current_slice_hits);
        auto [centroid_wire_v, centroid_drift_v] = this->calculateChargeWeightedCentroid(e, common::TPC_VIEW_V, current_slice_hits);
        auto [centroid_wire_w, centroid_drift_w] = this->calculateChargeWeightedCentroid(e, common::TPC_VIEW_W, current_slice_hits);

        std::vector<ImageProperties> view_img_props;
        view_img_props.emplace_back(centroid_wire_u, centroid_drift_u, _image_height_pixels, _image_width_pixels, _drift_step_cm_per_tick, _wire_pitch_u_cm, geo::kU);
        view_img_props.emplace_back(centroid_wire_v, centroid_drift_v, _image_height_pixels, _image_width_pixels, _drift_step_cm_per_tick, _wire_pitch_v_cm, geo::kV);
        view_img_props.emplace_back(centroid_wire_w, centroid_drift_w, _image_height_pixels, _image_width_pixels, _drift_step_cm_per_tick, _wire_pitch_w_cm, geo::kW);

        double nu_vtx_xyz[3] = {0.0, 0.0, 0.0};
        bool has_neutrino_vtx = false;
        art::Handle<std::vector<recob::PFParticle>> pfp_handle;
        if (e.getByLabel(fPFPproducer, pfp_handle)) {
            art::FindManyP<recob::Vertex> find_vertices(pfp_handle, e, fPFPproducer);
            if (find_vertices.isValid()) {
                for (const auto& pfp_ptr : slice_pfp_v) {
                    if (pfp_ptr.isNull()) continue;
                    const recob::PFParticle& pfp = *pfp_ptr;
                    if (pfp.IsPrimary() && (std::abs(pfp.PdgCode()) == 12 || std::abs(pfp.PdgCode()) == 14 || std::abs(pfp.PdgCode()) == 16)) {
                        std::vector<art::Ptr<recob::Vertex>> pfp_vertex_assn = find_vertices.at(pfp_ptr.key());
                        if (!pfp_vertex_assn.empty() && pfp_vertex_assn[0].isAvailable()) {
                            pfp_vertex_assn[0]->XYZ(nu_vtx_xyz);
                            has_neutrino_vtx = true;
                            break;
                        }
                    }
                }
            }
        }

        img_is_vertex_in_image = false;
        if (has_neutrino_vtx) {
            img_is_vertex_in_image = true;
            for (size_t view_idx = 0; view_idx < 3; ++view_idx) {
                const auto& props = view_img_props[view_idx];
                geo::View_t current_view = props.view();
                double vtx_wire_coord_proj;
                if (current_view == geo::kW) {
                    vtx_wire_coord_proj = nu_vtx_xyz[2];
                } else if (current_view == geo::kU) {
                    vtx_wire_coord_proj = nu_vtx_xyz[2] * std::cos(SIXTY_DEGREES_RAD) - nu_vtx_xyz[1] * std::sin(SIXTY_DEGREES_RAD);
                } else if (current_view == geo::kV) {
                    vtx_wire_coord_proj = nu_vtx_xyz[2] * std::cos(-SIXTY_DEGREES_RAD) - nu_vtx_xyz[1] * std::sin(-SIXTY_DEGREES_RAD);
                } else {
                    continue;
                }
                double vtx_drift_coord = nu_vtx_xyz[0];
                size_t r = props.row(vtx_drift_coord);
                size_t c = props.col(vtx_wire_coord_proj);
                if (r == static_cast<size_t>(-1) || c == static_cast<size_t>(-1) || r >= props.height() || c >= props.width()) {
                    img_is_vertex_in_image = false;
                    break;
                }
            }
        }

        std::vector<Image<float>> temp_raw_adc_imgs;
        std::vector<Image<int>> temp_cheated_reco_imgs;
        std::vector<Image<int>> temp_truth_primary_imgs;

        this->constructPixelImages(e, view_img_props, temp_raw_adc_imgs, temp_cheated_reco_imgs, temp_truth_primary_imgs);

        img_nrows = view_img_props.empty() ? 0 : view_img_props[0].height();
        img_ncols = view_img_props.empty() ? 0 : view_img_props[0].width();


        for (size_t view_idx = 0; view_idx < 3; ++view_idx) {
            if (view_idx < temp_raw_adc_imgs.size()) img_raw_adc_images[view_idx] = temp_raw_adc_imgs[view_idx].data();
            else if (view_idx < img_raw_adc_images.size()) img_raw_adc_images[view_idx].clear();

            if (fProcessMC) {
                if (view_idx < temp_cheated_reco_imgs.size()) img_cheated_reco_label_images[view_idx] = temp_cheated_reco_imgs[view_idx].data();
                else if (view_idx < img_cheated_reco_label_images.size()) img_cheated_reco_label_images[view_idx].clear();

                if (view_idx < temp_truth_primary_imgs.size()) img_truth_primary_label_images[view_idx] = temp_truth_primary_imgs[view_idx].data();
                else if (view_idx < img_truth_primary_label_images.size()) img_truth_primary_label_images[view_idx].clear();
            }
        }
        if (_image_tree_ptr) _image_tree_ptr->Fill();
    }

    std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        std::vector<art::Ptr<recob::Hit>> slice_hits;
        art::Handle<std::vector<recob::Cluster>> cluster_handle;
        if (!e.getByLabel(fCLSproducer, cluster_handle)) return slice_hits;
        art::FindManyP<recob::Hit> find_hits_from_clusters(cluster_handle, e, fCLSproducer);
        if (!find_hits_from_clusters.isValid()) return slice_hits;

        art::Handle<std::vector<recob::PFParticle>> pfp_handle;
        if (!e.getByLabel(fPFPproducer, pfp_handle)) return slice_hits;
        art::FindManyP<recob::Cluster> find_clusters_from_pfp(pfp_handle, e, fPFPproducer);
        if (!find_clusters_from_pfp.isValid()) return slice_hits;

        for (const auto& pfp_ptr : pfp_pxy_v) {
            if (pfp_ptr.isNull()) continue;
            if (pfp_ptr->IsPrimary()) continue;

            std::vector<art::Ptr<recob::Cluster>> associated_clusters = find_clusters_from_pfp.at(pfp_ptr.key());
            for (const auto& clus_ptr : associated_clusters) {
                if (clus_ptr.isNull()) continue;
                std::vector<art::Ptr<recob::Hit>> clus_hits_vec = find_hits_from_clusters.at(clus_ptr.key());
                slice_hits.insert(slice_hits.end(), clus_hits_vec.begin(), clus_hits_vec.end());
            }
        }
        std::sort(slice_hits.begin(), slice_hits.end());
        slice_hits.erase(std::unique(slice_hits.begin(), slice_hits.end()), slice_hits.end());
        return slice_hits;
    }

    std::pair<double, double> ImageAnalysis::calculateChargeWeightedCentroid(const art::Event& e, common::PandoraView target_pandora_view, const std::vector<art::Ptr<recob::Hit>>& hits) {
        double sum_charge = 0.0;
        double sum_weighted_wire_coord = 0.0;
        double sum_weighted_drift_coord = 0.0;
        if (hits.empty()) return {0.0, 0.0};

        geo::View_t target_geo_view = common::ConvertPandoraViewToGeoView(target_pandora_view);

        for (const auto& hit_ptr : hits) {
            if (hit_ptr.isNull()) continue;

            geo::View_t hit_geo_view = geo::kUnknown;
            if (hit_ptr->WireID().isValid) {
                 hit_geo_view = _geo_serv->View(hit_ptr->WireID().planeID());
            } else {
                 try { hit_geo_view = hit_ptr->View(); }
                 catch (const cet::exception& ) { continue; }
            }

            if (hit_geo_view != target_geo_view) continue;

            double charge = hit_ptr->Integral();
            if (charge <= 0) continue;

            TVector3 hit_pos_drift_wire = common::GetPandoraHitPosition(e, hit_ptr, target_pandora_view);
            sum_charge += charge;
            sum_weighted_drift_coord += hit_pos_drift_wire.X() * charge;
            sum_weighted_wire_coord += hit_pos_drift_wire.Z() * charge;
        }
        if (sum_charge == 0.0) return {0.0, 0.0};
        return {sum_weighted_wire_coord / sum_charge, sum_weighted_drift_coord / sum_charge};
    }

    void ImageAnalysis::constructPixelImages(const art::Event& e,
                                        const std::vector<ImageProperties>& view_img_props,
                                        std::vector<Image<float>>& out_raw_adc_imgs,
                                        std::vector<Image<int>>& out_cheated_reco_imgs,
                                        std::vector<Image<int>>& out_truth_primary_imgs
                                        ) {
        out_raw_adc_imgs.clear();
        if (fProcessMC) {
            out_cheated_reco_imgs.clear();
            out_truth_primary_imgs.clear();
        }

        for (const auto& props : view_img_props) {
            out_raw_adc_imgs.emplace_back(props);
            if (fProcessMC) {
                out_cheated_reco_imgs.emplace_back(props);
                out_truth_primary_imgs.emplace_back(props);
            }
        }

        auto wire_handle = e.getValidHandle<std::vector<recob::Wire>>(fWIREproducer);
        art::FindManyP<recob::Hit> wire_to_hits_assoc(wire_handle, e, fHITproducer);

        art::Handle<std::vector<simb::MCParticle>> mc_particle_handle;
        bool has_mc_particles = e.getByLabel(fMCPproducer, mc_particle_handle);
        
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> hit_to_mcparticle_data_assoc;
        std::vector<TruthPrimaryLabel> truth_primary_labels_for_event;
        std::vector<RecoLabel> reco_labels_for_event;
        std::map<int, size_t> trackid_to_mcparticle_idx;
        
        if (fProcessMC && has_mc_particles && mc_particle_handle.isValid() && !mc_particle_handle->empty()) {
            auto hit_handle = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
            hit_to_mcparticle_data_assoc = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hit_handle, e, fBackTrackerLabel);

            for (size_t i = 0; i < mc_particle_handle->size(); ++i) {
                trackid_to_mcparticle_idx[mc_particle_handle->at(i).TrackId()] = i;
            }
            truth_primary_labels_for_event = classifyTruthPrimaryLabels(e, fMCPproducer);
            reco_labels_for_event = classifyRecoLabels(e, fMCPproducer, fGammaThresholdGeV, fHadronThresholdGeV, fLeptonThresholdGeV);
        }
        
        for (size_t wire_idx = 0; wire_idx < wire_handle->size(); ++wire_idx) {
            const recob::Wire& wire = wire_handle->at(wire_idx);
            raw_channel_id_t channel = wire.Channel();
            std::vector<geo::WireID> wire_ids = _geo_serv->ChannelToWire(channel);
            if (wire_ids.empty()) continue;

            geo::WireID current_wid = wire_ids.front();
            geo::View_t current_view = current_wid.planeID().View;
            size_t view_idx = static_cast<size_t>(current_view);
            if (view_idx >= view_img_props.size()) continue;

            const ImageProperties& current_view_image_props = view_img_props[view_idx];
            const geo::WireGeo* wire_geometry = _geo_serv->WirePtr(current_wid);
            TVector3 wire_center_global_coords = wire_geometry->GetCenter();
            double wire_coord_projected;
            if (current_view == geo::kW) {
                wire_coord_projected = wire_center_global_coords.Z();
            } else if (current_view == geo::kU) {
                wire_coord_projected = wire_center_global_coords.Z() * std::cos(SIXTY_DEGREES_RAD) - wire_center_global_coords.Y() * std::sin(SIXTY_DEGREES_RAD);
            } else if (current_view == geo::kV) {
                wire_coord_projected = wire_center_global_coords.Z() * std::cos(-SIXTY_DEGREES_RAD) - wire_center_global_coords.Y() * std::sin(-SIXTY_DEGREES_RAD);
            } else {
                continue;
            }

            size_t image_col_idx = current_view_image_props.col(wire_coord_projected);
            if (image_col_idx == static_cast<size_t>(-1)) continue;

            std::vector<art::Ptr<recob::Hit>> hits_on_this_wire;
            if (wire_to_hits_assoc.isValid()) {
                hits_on_this_wire = wire_to_hits_assoc.at(wire_idx);
                std::sort(hits_on_this_wire.begin(), hits_on_this_wire.end(),
                          [](const art::Ptr<recob::Hit>& a, const art::Ptr<recob::Hit>& b) {
                              return a->StartTick() < b->StartTick();
                          });
            }
            size_t current_hit_search_iter = 0;

            for (const auto& roi_range : wire.SignalROI().get_ranges()) {
                const auto& adc_signal_values = roi_range.data();
                start_tick_roi = roi_range.begin_index();
                for (size_t adc_idx_in_roi = 0; adc_idx_in_roi < adc_signal_values.size(); ++adc_idx_in_roi) {
                    float adc_val = adc_signal_values[adc_idx_in_roi];
                    if (adc_val <= _adc_image_threshold) continue;

                    current_global_tick = start_tick_roi + adc_idx_in_roi;
                    double drift_coord_cm = _det_prop_serv->ConvertTicksToX(current_global_tick, current_wid.Plane, current_wid.TPC, current_wid.Cryostat);
                    size_t image_row_idx = current_view_image_props.row(drift_coord_cm);
                    if (image_row_idx == static_cast<size_t>(-1)) continue;

                    out_raw_adc_imgs[view_idx].set(image_row_idx, image_col_idx, adc_val, true);
                    
                    art::Ptr<recob::Hit> relevant_hit_ptr;
                    while (current_hit_search_iter < hits_on_this_wire.size() && hits_on_this_wire[current_hit_search_iter]->EndTick() < current_global_tick) {
                        current_hit_search_iter++;
                    }
                    if (current_hit_search_iter < hits_on_this_wire.size() &&
                        hits_on_this_wire[current_hit_search_iter]->StartTick() <= current_global_tick &&
                        current_global_tick < hits_on_this_wire[current_hit_search_iter]->EndTick()) {
                        relevant_hit_ptr = hits_on_this_wire[current_hit_search_iter];
                    }

                    if (fProcessMC) {
                        RecoLabel final_reco_label_for_pixel = RecoLabel::empty;
                        TruthPrimaryLabel final_truth_primary_label_for_pixel = TruthPrimaryLabel::Empty;

                        if (has_mc_particles && mc_particle_handle.isValid() && !mc_particle_handle->empty() && hit_to_mcparticle_data_assoc.isValid()) {
                            final_reco_label_for_pixel = RecoLabel::cosmic;
                            final_truth_primary_label_for_pixel = TruthPrimaryLabel::Cosmic;

                            if (!relevant_hit_ptr.isNull()) {
                                auto mcp_associations = hit_to_mcparticle_data_assoc.at(relevant_hit_ptr.key());
                                auto mcp_data_associations = hit_to_mcparticle_data_assoc.data(relevant_hit_ptr.key());
                                int best_mcp_track_id = -1;

                                if (!mcp_associations.empty() && !mcp_data_associations.empty()) {
                                    for (size_t mcp_assoc_idx = 0; mcp_assoc_idx < mcp_associations.size(); ++mcp_assoc_idx) {
                                        if (mcp_data_associations[mcp_assoc_idx] && mcp_data_associations[mcp_assoc_idx]->isMaxIDE) {
                                            best_mcp_track_id = mcp_associations[mcp_assoc_idx]->TrackId();
                                            break;
                                        }
                                    }
                                }

                                if (best_mcp_track_id != -1) {
                                    auto it_mcp = trackid_to_mcparticle_idx.find(best_mcp_track_id);
                                    if (it_mcp != trackid_to_mcparticle_idx.end()) {
                                        size_t mc_part_idx_in_event_list = it_mcp->second;
                                        if (mc_part_idx_in_event_list < reco_labels_for_event.size()) {
                                            final_reco_label_for_pixel = reco_labels_for_event[mc_part_idx_in_event_list];
                                        }
                                        if (mc_part_idx_in_event_list < truth_primary_labels_for_event.size()) {
                                            final_truth_primary_label_for_pixel = truth_primary_labels_for_event[mc_part_idx_in_event_list];
                                        }
                                    }
                                }
                            }
                        }
                        out_cheated_reco_imgs[view_idx].set(image_row_idx, image_col_idx, static_cast<int>(final_reco_label_for_pixel), false);
                        out_truth_primary_imgs[view_idx].set(image_row_idx, image_col_idx, static_cast<int>(final_truth_primary_label_for_pixel), false);
                    }
                }
            }
        }
    }
    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}
