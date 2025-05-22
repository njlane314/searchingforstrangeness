#ifndef ANALYSIS_IMAGE_CXX
#define ANALYSIS_IMAGE_CXX

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
#include <optional>

#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TVector3.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "AnalysisToolBase.h"
#include "../CommonDefs/Pandora.h"
#include "../CommonDefs/Types.h"

namespace analysis 
{
    class ImageProperties {
    public:
        ImageProperties() = default;
        ImageProperties(double center_x, double center_y, size_t width, size_t height, double pixel_h, double pixel_w, geo::View_t view)
        : center_x_(center_x), center_y_(center_y), height_(height), width_(width), pixel_w_(pixel_w), pixel_h_(pixel_h), view_(view) {
            origin_x_ = center_x - (width * pixel_w) / 2.0;
            origin_y_ = center_y - (height * pixel_h) / 2.0;
        }
        size_t index(size_t row, size_t col) const {
            if (row >= height_ || col >= width_) return static_cast<size_t>(-1);
            return col * height_ + row;
        }
        size_t col(double x) const {
            if (x < origin_x_ || x >= origin_x_ + width_ * pixel_w_) return static_cast<size_t>(-1);
            return static_cast<size_t>((x - origin_x_) / pixel_w_);
        }
        size_t row(double y) const {
            if (y < origin_y_ || y >= origin_y_ + height_ * pixel_h_) return static_cast<size_t>(-1);
            return static_cast<size_t>((y - origin_y_) / pixel_h_);
        }
        size_t height() const { return height_; }
        size_t width() const { return width_; }
        double pixel_w() const { return pixel_w_; }
        double pixel_h() const { return pixel_h_; }
        geo::View_t view() const { return view_; }
        double origin_x() const { return origin_x_; }
        double origin_y() const { return origin_y_; }
        double max_x() const { return origin_x_ + width_ * pixel_w_; }
        double max_y() const { return origin_y_ + height_ * pixel_h_; }
    private:
        double center_x_, center_y_;
        double origin_x_, origin_y_;
        size_t height_, width_;
        double pixel_w_, pixel_h_;
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
        std::vector<T> data() const {
            return pixels_;
        }
        geo::View_t view() const { return prop_.view(); }
        size_t height() const { return prop_.height(); }
        size_t width() const { return prop_.width(); }
    private:
        ImageProperties prop_;
        std::vector<T> pixels_;
    };

    class TruthLabelClassifier {
    public:
        enum class TruthPrimaryLabel {
            Empty = 0,
            Cosmic,
            Electron,
            Muon,
            ChargedPion,
            NeutralPion,
            Proton,
            ChargedKaon,
            NeutralKaon,
            Lambda,
            ChargedSigma,
            Other
        };

        static inline const std::array<std::string, 12> truth_primary_label_names = {
            "Empty", "Cosmic", "Electron", "Muon", "ChargedPion", "NeutralPion", "Proton", "ChargedKaon",
            "NeutralKaon", "Lambda", "ChargedSigma", "Other"
        };

        explicit TruthLabelClassifier(const art::InputTag& MCPproducer)
        : fMCPproducer(MCPproducer) {}

        TruthPrimaryLabel getTruthPrimaryLabelFromPDG(int pdg) const {
            pdg = std::abs(pdg);
            if (pdg == 11) return TruthPrimaryLabel::Electron;
            if (pdg == 13) return TruthPrimaryLabel::Muon;
            if (pdg == 211) return TruthPrimaryLabel::ChargedPion;
            if (pdg == 111) return TruthPrimaryLabel::NeutralPion;
            if (pdg == 2212) return TruthPrimaryLabel::Proton;
            if (pdg == 321) return TruthPrimaryLabel::ChargedKaon;
            if (pdg == 311 || pdg == 130 || pdg == 310) return TruthPrimaryLabel::NeutralKaon;
            if (pdg == 3122) return TruthPrimaryLabel::Lambda;
            if (pdg == 3222 || pdg == 3112 || pdg == 3212) return TruthPrimaryLabel::ChargedSigma;
            return TruthPrimaryLabel::Other;
        }

        void assignLabelToProgenyRecursively(
            size_t particle_index,
            const std::vector<simb::MCParticle>& particles,
            std::vector<TruthPrimaryLabel>& particle_labels,
            const std::unordered_map<int, size_t>& track_id_to_index,
            TruthPrimaryLabel primary_label_to_assign
        ) const {
            if (particle_index >= particles.size() || particle_index >= particle_labels.size()) {
                return;
            }
            particle_labels[particle_index] = primary_label_to_assign;
            const auto& particle = particles[particle_index];

            for (int daughter_idx = 0; daughter_idx < particle.NumberDaughters(); ++daughter_idx) {
                int daughter_track_id = particle.Daughter(daughter_idx);
                auto it = track_id_to_index.find(daughter_track_id);
                if (it != track_id_to_index.end()) {
                    if (it->second < particles.size()) {
                        assignLabelToProgenyRecursively(it->second, particles, particle_labels, track_id_to_index, primary_label_to_assign);
                    }
                }
            }
        }

        std::vector<TruthPrimaryLabel> classifyParticles(
            const art::Event& event
        ) const {
            const auto particle_collection_handle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
            const auto& particles = *particle_collection_handle;

            std::unordered_map<int, size_t> track_id_to_vector_index;
            for (size_t i = 0; i < particles.size(); ++i) {
                track_id_to_vector_index[particles[i].TrackId()] = i;
            }

            std::vector<TruthPrimaryLabel> classified_particle_labels(particles.size(), TruthPrimaryLabel::Empty);

            for (size_t i = 0; i < particles.size(); ++i) {
                if (particles[i].Mother() == 0) {
                    if (auto it = track_id_to_vector_index.find(particles[i].TrackId()); it != track_id_to_vector_index.end()) {
                        TruthPrimaryLabel initial_label = getTruthPrimaryLabelFromPDG(particles[i].PdgCode());
                        assignLabelToProgenyRecursively(it->second, particles, classified_particle_labels, track_id_to_vector_index, initial_label);
                    }
                }
            }
            return classified_particle_labels;
        }

    private:
        art::InputTag fMCPproducer;
    };

    class RecoLabelClassifier {
    public:
        enum class ReconstructionLabel {
            Empty,
            Cosmic,
            MIP,
            HIP,
            Shower,
            Michel,
            Diffuse,
            Invisible
        };

        static inline const std::array<std::string, 8> reco_label_names = {
            "Empty", "Cosmic", "MIP", "HIP", "Shower", "Michel", "Diffuse", "Invisible"
        };

        explicit RecoLabelClassifier(
            const art::InputTag& mc_producer_tag,
            double gamma_threshold,
            double hadron_threshold,
            bool use_cheat)
        : fMCPproducer(mc_producer_tag),
        fGammaThreshold(gamma_threshold),
        fHadronThreshold(hadron_threshold),
        fUseCheat(use_cheat) {}

        std::pair<ReconstructionLabel, ReconstructionLabel> getReconstructionLabelAndPropagation(
            const std::vector<simb::MCParticle>& particles,
            const simb::MCParticle& particle_to_label,
            ReconstructionLabel label_from_parent,
            const std::map<int, size_t>& track_id_to_index
        ) const {
            if (label_from_parent != ReconstructionLabel::Empty) {
                return {label_from_parent, label_from_parent};
            }

            ReconstructionLabel determined_label = ReconstructionLabel::Invisible;
            ReconstructionLabel propagated_label_for_children = ReconstructionLabel::Empty;

            int pdg_code = particle_to_label.PdgCode();
            double momentum = particle_to_label.P();
            const std::string& start_process = particle_to_label.Process();
            const std::string& end_process = particle_to_label.EndProcess();
            int parent_track_id = particle_to_label.Mother();
            int parent_pdg_code = 0;

            if (parent_track_id != 0) {
                auto it = track_id_to_index.find(parent_track_id);
                if (it != track_id_to_index.end()) {
                    if (it->second < particles.size()) {
                        parent_pdg_code = particles[it->second].PdgCode();
                    }
                }
            }

            if (std::abs(pdg_code) == 211 || std::abs(pdg_code) == 13) {
                determined_label = ReconstructionLabel::MIP;
            } else if (std::abs(pdg_code) == 321 || (std::abs(pdg_code) == 2212 && momentum >= fHadronThreshold)) {
                determined_label = ReconstructionLabel::HIP;
            } else if (std::abs(pdg_code) == 11) {
                if (start_process == "primary") {
                    determined_label = ReconstructionLabel::Shower;
                    propagated_label_for_children = ReconstructionLabel::Shower;
                } else if (std::abs(parent_pdg_code) == 13 &&
                        (start_process == "muMinusCaptureAtRest" ||
                            start_process == "muPlusCaptureAtRest" ||
                            start_process == "Decay")) {
                    determined_label = ReconstructionLabel::Michel;
                    propagated_label_for_children = ReconstructionLabel::Michel;
                } else if (start_process == "conv" || end_process == "conv" ||
                        start_process == "compt" || end_process == "compt") {
                    if (momentum >= fGammaThreshold) {
                        determined_label = ReconstructionLabel::Shower;
                        propagated_label_for_children = ReconstructionLabel::Shower;
                    } else {
                        determined_label = ReconstructionLabel::Diffuse;
                    }
                } else {
                    determined_label = ReconstructionLabel::Diffuse;
                }
            } else if (pdg_code == 22) {
                if (start_process == "conv" || end_process == "conv" ||
                    start_process == "compt" || end_process == "compt" ||
                    start_process == "primary") {
                    if (momentum >= fGammaThreshold) {
                        determined_label = ReconstructionLabel::Shower;
                        propagated_label_for_children = ReconstructionLabel::Shower;
                    } else {
                        determined_label = ReconstructionLabel::Diffuse;
                    }
                } else {
                    determined_label = ReconstructionLabel::Diffuse;
                }
            } else if (std::abs(pdg_code) == 2212 && momentum < fHadronThreshold) {
                determined_label = ReconstructionLabel::Diffuse;
            }
            return {determined_label, propagated_label_for_children};
        }

        void assignLabelToProgenyRecursively(
            size_t particle_index,
            const std::vector<simb::MCParticle>& particles,
            ReconstructionLabel label_from_parent,
            std::vector<ReconstructionLabel>& particle_labels,
            const std::map<int, size_t>& track_id_to_index
        ) const {
            if (particle_index >= particles.size() || particle_index >= particle_labels.size()) {
                return;
            }
            const auto& current_particle = particles[particle_index];
            auto [label_for_current, label_for_children] = getReconstructionLabelAndPropagation(
                particles, current_particle, label_from_parent, track_id_to_index
            );
            particle_labels[particle_index] = label_for_current;

            for (int i = 0; i < current_particle.NumberDaughters(); ++i) {
                int daughter_track_id = current_particle.Daughter(i);
                auto it = track_id_to_index.find(daughter_track_id);
                if (it != track_id_to_index.end()) {
                    if (it->second < particles.size()) {
                        assignLabelToProgenyRecursively(it->second, particles, label_for_children, particle_labels, track_id_to_index);
                    }
                }
            }
        }

        std::vector<ReconstructionLabel> classifyParticlesFromReco(const art::Event& event) const {
            return {};
        }

        std::vector<ReconstructionLabel> classifyParticles(const art::Event& event) const {
            if (fUseCheat) {
                const auto particle_collection_handle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
                const auto& particles = *particle_collection_handle;

                std::map<int, size_t> track_id_to_vector_index;
                for (size_t i = 0; i < particles.size(); ++i) {
                    track_id_to_vector_index[particles[i].TrackId()] = i;
                }

                std::vector<ReconstructionLabel> classified_particle_labels(particles.size(), ReconstructionLabel::Empty);

                for (size_t i = 0; i < particles.size(); ++i) {
                    if (particles[i].Mother() == 0) {
                        assignLabelToProgenyRecursively(i, particles, ReconstructionLabel::Empty, classified_particle_labels, track_id_to_vector_index);
                    }
                }
                return classified_particle_labels;
            } else {
                return classifyParticlesFromReco(event);
            }
        }

        art::InputTag fMCPproducer;
        double fGammaThreshold;
        double fHadronThreshold;
        bool fUseCheat;
    };


    class ImageAnalysis : public AnalysisToolBase {
    public:
        explicit ImageAnalysis(fhicl::ParameterSet const& p);
        virtual ~ImageAnalysis() = default;

        void configure(const fhicl::ParameterSet& p) override;
        void analyseEvent(art::Event const& e, bool _is_data) override {}
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) override;
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;

    private:
        art::InputTag fPFPproducer;
        art::InputTag fCLSproducer;
        art::InputTag fHITproducer;
        art::InputTag fWIREproducer;
        art::InputTag fMCPproducer;
        std::string fBackTrackerLabel;

        bool fProcessMC;
        bool fUseCheatRecoLabels;

        double fGammaThreshold;
        double fHadronThreshold;
        double fLeptonThreshold;

        int _image_width;
        int _image_height;
        float _adc_image_threshold;

        const geo::GeometryCore* _geo;
        const detinfo::DetectorProperties* _detp;

        float _drift_step;
        float _wire_pitch_u;
        float _wire_pitch_v;
        float _wire_pitch_w;

        TTree* _tree = nullptr;

        std::vector<float> _raw_image_u;
        std::vector<float> _raw_image_v;
        std::vector<float> _raw_image_w;
        std::vector<int> _reco_image_u;
        std::vector<int> _reco_image_v;
        std::vector<int> _reco_image_w;
        std::vector<int> _true_image_u;
        std::vector<int> _true_image_v;
        std::vector<int> _true_image_w;

        std::unique_ptr<RecoLabelClassifier> fRecoClassifier;
        std::unique_ptr<TruthLabelClassifier> fTruthClassifier;

        std::vector<art::Ptr<recob::Hit>> collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);
        std::pair<double, double> calculateChargeCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits);
        void constructPixelImages(const art::Event& e, const std::vector<ImageProperties>& properties,
                                std::vector<Image<float>>& raw_images,
                                std::vector<Image<int>>& reco_images,
                                std::vector<Image<int>>& true_images);
    };

    ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void ImageAnalysis::configure(const fhicl::ParameterSet& p) {
        fPFPproducer = p.get<art::InputTag>("PFPproducer");
        fCLSproducer = p.get<art::InputTag>("CLSproducer");
        fHITproducer = p.get<art::InputTag>("HITproducer");
        fWIREproducer = p.get<art::InputTag>("WIREproducer");
        fMCPproducer = p.get<art::InputTag>("MCPproducer");
        fBackTrackerLabel = p.get<std::string>("BackTrackerLabel", "gaushit");

        fGammaThreshold = p.get<double>("GammaThreshold", 0.1);
        fHadronThreshold = p.get<double>("HadronThreshold", 0.1);
        fLeptonThreshold = p.get<double>("LeptonThreshold", 0.1);

        fProcessMC = p.get<bool>("ProcessMC", true);
        fUseCheatRecoLabels = p.get<bool>("UseCheatRecoLabels", true);

        _image_width = p.get<int>("ImageWidth", 512);
        _image_height = p.get<int>("ImageHeight", 512);
        _adc_image_threshold = p.get<float>("ADCthreshold", 4.0);

        _geo = art::ServiceHandle<geo::Geometry>()->provider();
        _detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();

        double tick_period = clock->TPCClock().TickPeriod();
        double drift_velocity = _detp->DriftVelocity();
        _drift_step = tick_period * drift_velocity * 1e1;
        _wire_pitch_u = _geo->WirePitch(geo::kU);
        _wire_pitch_v = _geo->WirePitch(geo::kV);
        _wire_pitch_w = _geo->WirePitch(geo::kW);

        if (fProcessMC) {
            fTruthClassifier = std::make_unique<TruthLabelClassifier>(fMCPproducer);
            fRecoClassifier = std::make_unique<RecoLabelClassifier>(fMCPproducer, fGammaThreshold, fHadronThreshold, fUseCheatRecoLabels);
        }
    }

    void ImageAnalysis::setBranches(TTree* _tree_ptr) {
        _tree = _tree_ptr;
        _tree->Branch("raw_image_u", &_raw_image_u);
        _tree->Branch("raw_image_v", &_raw_image_v);
        _tree->Branch("raw_image_w", &_raw_image_w);
        _tree->Branch("reco_image_u", &_reco_image_u);
        _tree->Branch("reco_image_v", &_reco_image_v);
        _tree->Branch("reco_image_w", &_reco_image_w);
        _tree->Branch("true_image_u", &_true_image_u);
        _tree->Branch("true_image_v", &_true_image_v);
        _tree->Branch("true_image_w", &_true_image_w);
    }

    void ImageAnalysis::resetTTree(TTree* _tree_ptr) {
        _raw_image_u.clear();
        _raw_image_v.clear();
        _raw_image_w.clear();
        _reco_image_u.clear();
        _reco_image_v.clear();
        _reco_image_w.clear();
        _true_image_u.clear();
        _true_image_v.clear();
        _true_image_w.clear();
    }

    void ImageAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits = this->collectSliceHits(e, slice_pfp_v);
        auto [centroid_wire_u, centroid_drift_u] = this->calculateChargeCentroid(e, common::TPC_VIEW_U, neutrino_hits);
        auto [centroid_wire_v, centroid_drift_v] = this->calculateChargeCentroid(e, common::TPC_VIEW_V, neutrino_hits);
        auto [centroid_wire_w, centroid_drift_w] = this->calculateChargeCentroid(e, common::TPC_VIEW_W, neutrino_hits);
        std::vector<ImageProperties> properties;
        properties.emplace_back(centroid_wire_u, centroid_drift_u, _image_height, _image_width, _drift_step, _wire_pitch_u, geo::kU);
        properties.emplace_back(centroid_wire_v, centroid_drift_v, _image_height, _image_width, _drift_step, _wire_pitch_v, geo::kV);
        properties.emplace_back(centroid_wire_w, centroid_drift_w, _image_height, _image_width, _drift_step, _wire_pitch_w, geo::kW);

        std::vector<Image<float>> out_raw_images;
        std::vector<Image<int>> out_reco_images;
        std::vector<Image<int>> out_true_images;
        this->constructPixelImages(e, properties, out_raw_images, out_reco_images, out_true_images);

        _raw_image_u.clear();
        _raw_image_v.clear();
        _raw_image_w.clear();
        _reco_image_u.clear();
        _reco_image_v.clear();
        _reco_image_w.clear();
        _true_image_u.clear();
        _true_image_v.clear();
        _true_image_w.clear();

        _raw_image_u = out_raw_images[0].data();  // U view (geo::kU)
        _raw_image_v = out_raw_images[1].data();  // V view (geo::kV)
        _raw_image_w = out_raw_images[2].data();  // W view (geo::kW)
        _reco_image_u = out_reco_images[0].data();
        _reco_image_v = out_reco_images[1].data();
        _reco_image_w = out_reco_images[2].data();
        _true_image_u = out_true_images[0].data();
        _true_image_v = out_true_images[1].data();
        _true_image_w = out_true_images[2].data();

        if (_tree) _tree->Fill();
    }

    std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits;
        auto clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
        for (const auto& pfp : pfp_pxy_v) {
            if (pfp->IsPrimary()) continue;
            for (auto ass_clus : pfp.get<recob::Cluster>()) {
                auto clus_hit_v = clus_proxy[ass_clus.key()].get<recob::Hit>();
                neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
            }
        }
        return neutrino_hits;
    }

    std::pair<double, double> ImageAnalysis::calculateChargeCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits) {
        double sum_charge = 0.0;
        double sum_wire = 0.0;
        double sum_drift = 0.0;
        for (const auto& hit : hits) {
            if (common::GetPandoraView(hit) != view) continue;
            double charge = hit->Integral();
            TVector3 hit_pos = common::GetPandoraHitPosition(e, hit, view);
            sum_charge += charge;
            sum_wire += hit_pos.Z() * charge;
            sum_drift += hit_pos.X() * charge;
        }
        if (sum_charge == 0.0) return {0.0, 0.0};
        return {sum_wire / sum_charge, sum_drift / sum_charge};
    }

    void ImageAnalysis::constructPixelImages(const art::Event& e,
                                            const std::vector<ImageProperties>& properties,
                                            std::vector<Image<float>>& out_raw_images,
                                            std::vector<Image<int>>& out_reco_images,
                                            std::vector<Image<int>>& out_true_images) {
        out_raw_images.clear();
        out_reco_images.clear();
        out_true_images.clear();

        for (const auto& prop : properties) {
            Image<float> raw_image(prop);
            raw_image.clear(0.0);
            out_raw_images.push_back(std::move(raw_image));

            Image<int> reco_image(prop);
            reco_image.clear(static_cast<int>(RecoLabelClassifier::ReconstructionLabel::Empty));
            out_reco_images.push_back(std::move(reco_image));

            Image<int> true_image(prop);
            true_image.clear(static_cast<int>(TruthLabelClassifier::TruthPrimaryLabel::Empty));
            out_true_images.push_back(std::move(true_image));
        }

        auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(fWIREproducer);
        auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);

        art::Handle<std::vector<simb::MCParticle>> mcpHandle;
        bool hasMCParticles = e.getByLabel(fMCPproducer, mcpHandle);

        art::FindManyP<recob::Hit> wire_hit_assoc(wireHandle, e, fHITproducer);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hitHandle, e, fBackTrackerLabel);

        std::vector<TruthLabelClassifier::TruthPrimaryLabel> classified_true_labels;
        std::vector<RecoLabelClassifier::ReconstructionLabel> classified_reco_labels;

        if (fProcessMC && hasMCParticles && mcpHandle.isValid() && fTruthClassifier && fRecoClassifier) {
            classified_true_labels = fTruthClassifier->classifyParticles(e);
            classified_reco_labels = fRecoClassifier->classifyParticles(e);
        }

        std::map<int, size_t> trackid_to_index;
        if (fProcessMC && hasMCParticles && mcpHandle.isValid()) {
            for (size_t i = 0; i < mcpHandle->size(); ++i) {
                trackid_to_index[mcpHandle->at(i).TrackId()] = i;
            }
        }

        for (size_t wire_idx = 0; wire_idx < wireHandle->size(); ++wire_idx) {
            const auto& wire = wireHandle->at(wire_idx);
            auto ch_id = wire.Channel();
            std::vector<geo::WireID> wire_ids = _geo->ChannelToWire(ch_id);
            if (wire_ids.empty()) continue;
            geo::View_t view = _geo->View(wire_ids.front().planeID());
            size_t view_idx = static_cast<size_t>(view); 

            const geo::WireGeo* wire_geo = _geo->WirePtr(wire_ids.front());
            TVector3 center = wire_geo->GetCenter();
            TVector3 wire_center(center.X(), center.Y(), center.Z());
            double wire_coord = (view == geo::kW) ? wire_center.Z() :
                                (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));

            auto hits_for_wire = wire_hit_assoc.at(wire_idx);
            std::vector<art::Ptr<recob::Hit>> sorted_hits = hits_for_wire;


            std::sort(sorted_hits.begin(), sorted_hits.end(), [](const art::Ptr<recob::Hit>& a, const art::Ptr<recob::Hit>& b) {
                return a->StartTick() < b->StartTick();
            });

            size_t current_hit_idx_in_sorted_list = 0;

            for (const auto& range : wire.SignalROI().get_ranges()) {
                const auto& adcs = range.data();
                int start_tick = range.begin_index();
                for (size_t idx_in_adc_range = 0; idx_in_adc_range < adcs.size(); ++idx_in_adc_range) {
                    int tick = start_tick + idx_in_adc_range;
                    double x_drift_pos = _detp->ConvertTicksToX(static_cast<double>(tick), wire_ids.front().planeID());
                    size_t row = properties[view_idx].row(x_drift_pos);
                    size_t col = properties[view_idx].col(wire_coord);
                    if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                    RecoLabelClassifier::ReconstructionLabel current_reco_label_for_pixel;
                    TruthLabelClassifier::TruthPrimaryLabel current_true_label_for_pixel;

                    if (fProcessMC && hasMCParticles && mcpHandle.isValid() && fTruthClassifier && fRecoClassifier) {
                        current_reco_label_for_pixel = RecoLabelClassifier::ReconstructionLabel::Invisible;
                        current_true_label_for_pixel = TruthLabelClassifier::TruthPrimaryLabel::Other;

                        while (current_hit_idx_in_sorted_list < sorted_hits.size() && sorted_hits[current_hit_idx_in_sorted_list]->EndTick() <= tick) {
                            ++current_hit_idx_in_sorted_list;
                        }

                        if (current_hit_idx_in_sorted_list < sorted_hits.size() &&
                            sorted_hits[current_hit_idx_in_sorted_list]->StartTick() <= tick &&
                            tick < sorted_hits[current_hit_idx_in_sorted_list]->EndTick()) {
                            const art::Ptr<recob::Hit>& matched_hit = sorted_hits[current_hit_idx_in_sorted_list];
                            
                            std::vector<art::Ptr<simb::MCParticle>> mcp_particles_ass_to_hit;
                            std::vector<anab::BackTrackerHitMatchingData const*> bkth_data_ass_to_hit;
                            mcp_bkth_assoc.get(matched_hit.key(), mcp_particles_ass_to_hit, bkth_data_ass_to_hit);


                            if (!bkth_data_ass_to_hit.empty()) {
                                for (size_t i_bkth = 0; i_bkth < bkth_data_ass_to_hit.size(); ++i_bkth) {
                                    if (bkth_data_ass_to_hit[i_bkth] && bkth_data_ass_to_hit[i_bkth]->isMaxIDE) {
                                        int track_id = mcp_particles_ass_to_hit[i_bkth]->TrackId();
                                        auto it_trackid = trackid_to_index.find(track_id);
                                        if (it_trackid != trackid_to_index.end()) {
                                            size_t particle_mcp_idx = it_trackid->second;
                                            if (particle_mcp_idx < classified_reco_labels.size()) {
                                                current_reco_label_for_pixel = classified_reco_labels[particle_mcp_idx];
                                            }
                                            if (particle_mcp_idx < classified_true_labels.size()) {
                                                current_true_label_for_pixel = classified_true_labels[particle_mcp_idx];
                                            }
                                        }
                                        break; 
                                    }
                                }
                            } else {
                                current_reco_label_for_pixel = RecoLabelClassifier::ReconstructionLabel::Cosmic;
                                current_true_label_for_pixel = TruthLabelClassifier::TruthPrimaryLabel::Cosmic;
                            }
                        }
                    } else {
                        current_reco_label_for_pixel = RecoLabelClassifier::ReconstructionLabel::Cosmic;
                        current_true_label_for_pixel = TruthLabelClassifier::TruthPrimaryLabel::Cosmic;
                    }

                    if (adcs[idx_in_adc_range] > _adc_image_threshold) {
                        out_raw_images[view_idx].set(row, col, adcs[idx_in_adc_range]);
                        out_reco_images[view_idx].set(row, col, static_cast<int>(current_reco_label_for_pixel), false);
                        out_true_images[view_idx].set(row, col, static_cast<int>(current_true_label_for_pixel), false);
                    }
                }
            }
        }
    }

    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif