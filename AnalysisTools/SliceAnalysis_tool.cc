#ifndef ANALYSIS_SLICE_CXX
#define ANALYSIS_SLICE_CXX

#include <iostream>
#include <map>
#include <unordered_set>
#include "AnalysisToolBase.h"
#include "larcore/Geometry/Geometry.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "TVector3.h"
#include "../CommonFunctions/Backtracking.h"
#include "../CommonFunctions/Scores.h"
#include "../CommonFunctions/Corrections.h"
#include "../CommonFunctions/Scatters.h"
#include "../CommonFunctions/Pandora.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "../EventClassifier.h"
#include <Eigen/Dense>
#include "lardataobj/MCBase/MCShower.h"

namespace analysis 
{
    class SliceAnalysis : public AnalysisToolBase {
    public:
        explicit SliceAnalysis(const fhicl::ParameterSet& pset);
        ~SliceAnalysis() = default;

        SliceAnalysis(const SliceAnalysis&) = delete;
        SliceAnalysis(SliceAnalysis&&) = delete;
        SliceAnalysis& operator=(const SliceAnalysis&) = delete;
        SliceAnalysis& operator=(SliceAnalysis&&) = delete;

        void configure(const fhicl::ParameterSet& pset) override;
        void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) override;
        void analyseEvent(const art::Event& event, bool is_data) override {}
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        art::InputTag _PandoraModuleLabel;
        art::InputTag _HitModuleLabel;
        art::InputTag _BacktrackModuleLabel;
        art::InputTag _PFParticleModuleLabel;
        art::InputTag _SpacePointModuleLabel;
        art::InputTag _MCPproducer;
        art::InputTag _CLSproducer;

        struct Slice {
            float hit_completeness = std::numeric_limits<float>::lowest();
            float hit_purity = std::numeric_limits<float>::lowest();
            float charge_completeness = std::numeric_limits<float>::lowest();
            float charge_purity = std::numeric_limits<float>::lowest();
            float pattern_hit_purity = std::numeric_limits<float>::lowest();
            float pattern_hit_completeness = std::numeric_limits<float>::lowest();
            float pattern_charge_purity = std::numeric_limits<float>::lowest();
            float pattern_charge_completeness = std::numeric_limits<float>::lowest();
            int identifier = -1;
            int num_hits = -1;
            float center_x = std::numeric_limits<float>::lowest();
            float center_y = std::numeric_limits<float>::lowest();
            float center_z = std::numeric_limits<float>::lowest();
            float charge = std::numeric_limits<float>::lowest();
            float total_charge = std::numeric_limits<float>::lowest();
            float pca_eigenvalues[3] = {std::numeric_limits<float>::lowest(), 
                            std::numeric_limits<float>::lowest(), 
                            std::numeric_limits<float>::lowest()};
            float pca_principal_dir[3] = {0.0f, 0.0f, 0.0f};
        };

        Slice _event_slice;
        int _total_true_hits;
        std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> _assoc_mc_part;
        std::unique_ptr<signature::EventClassifier> _eventClassifier;
        int _total_neutrino_hits;
        int _total_pattern_hits;
        double _total_neutrino_charge;
        double _total_pattern_charge;
        std::unordered_set<int> _neutrino_track_ids;
        std::unordered_set<int> _pattern_track_ids;

        std::vector<art::Ptr<recob::Hit>> getNeutrinoHits(
            const art::Event& e,
            const std::vector<common::ProxyPfpElem_t>& neutrino_slice,
            const art::InputTag& clusterProducerTag);
        void fillSlice(const art::Event& event, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v);
    };

    SliceAnalysis::SliceAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void SliceAnalysis::configure(const fhicl::ParameterSet& pset) {
        _PandoraModuleLabel = pset.get<art::InputTag>("PandoraModuleLabel", "pandora");
        _HitModuleLabel = pset.get<art::InputTag>("HitModuleLabel", "gaushit");
        _BacktrackModuleLabel = pset.get<art::InputTag>("BacktrackModuleLabel", "gaushitTruthMatch");
        _PFParticleModuleLabel = pset.get<art::InputTag>("PFParticleModuleLabel", "pandora");
        _SpacePointModuleLabel = pset.get<art::InputTag>("SpacePointModuleLabel", "pandora");
        _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant"); 
        _CLSproducer = pset.get<art::InputTag>("CLSproducer", "pandora");
        _eventClassifier = std::make_unique<signature::EventClassifier>(pset.get<fhicl::ParameterSet>("EventClassifier"));
    }

    void SliceAnalysis::analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) {
        if (is_data) return;

        auto mcp_h = event.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
        _neutrino_track_ids.clear();
        for (const auto& mcp : *mcp_h) {
            _neutrino_track_ids.insert(mcp.TrackId());
        }

        const auto& pattern = _eventClassifier->getPattern(event);

        _pattern_track_ids.clear();
        for (const auto& [type, sig] : pattern) {
            for (const auto& particle : sig) {
                _pattern_track_ids.insert(particle->TrackId());
            }
        }

        auto const& hit_handle = event.getValidHandle<std::vector<recob::Hit>>(_HitModuleLabel);
        _assoc_mc_part = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_handle, event, _BacktrackModuleLabel);

        _total_neutrino_hits = 0;
        _total_pattern_hits = 0;
        _total_neutrino_charge = 0.0;
        _total_pattern_charge = 0.0;

        for (size_t i = 0; i < hit_handle->size(); ++i) {
            const auto& hit = (*hit_handle)[i];
            const auto& matched_mc_part = _assoc_mc_part->at(i);
            const auto& matched_data = _assoc_mc_part->data(i);
            double hit_charge = hit.Integral();
            double neutrino_charge = 0.0;
            double pattern_charge = 0.0;
            bool is_neutrino_hit = false;
            bool is_pattern_hit = false;
            for (size_t j = 0; j < matched_mc_part.size(); ++j) {
                int track_id = matched_mc_part[j]->TrackId();
                if (_neutrino_track_ids.count(track_id)) {
                    neutrino_charge += hit_charge * matched_data[j]->ideFraction;
                    if (matched_data[j]->isMaxIDE == 1) {
                        is_neutrino_hit = true;
                    }
                }
                if (_pattern_track_ids.count(track_id)) {
                    pattern_charge += hit_charge * matched_data[j]->ideFraction;
                    if (matched_data[j]->isMaxIDE == 1) {
                        is_pattern_hit = true;
                    }
                }
            }
            _total_neutrino_charge += neutrino_charge;
            _total_pattern_charge += pattern_charge;
            if (is_neutrino_hit) {
                ++_total_neutrino_hits;
            }
            if (is_pattern_hit) {
                ++_total_pattern_hits;
            }
        }

        this->fillSlice(event, slice_pfp_v);
    }

    std::vector<art::Ptr<recob::Hit>> SliceAnalysis::getNeutrinoHits(
        const art::Event& e,
        const std::vector<common::ProxyPfpElem_t>& neutrino_slice,
        const art::InputTag& clusterProducerTag) {

        std::vector<art::Ptr<recob::Hit>> neutrino_hits;
        auto clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(
            e, clusterProducerTag, proxy::withAssociated<recob::Hit>(clusterProducerTag));
        for (const auto& pfp : neutrino_slice) {
            if (pfp->IsPrimary()) continue;
            for (const auto& ass_clus : pfp.get<recob::Cluster>()) {
                auto clus_hit_v = clus_proxy[ass_clus.key()].get<recob::Hit>();
                neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
            }
        }
        return neutrino_hits;
    }

    void SliceAnalysis::fillSlice(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v) {
        std::vector<art::Ptr<recob::Hit>> slice_hits = getNeutrinoHits(e, slice_pfp_v, _CLSproducer);

        _event_slice.num_hits = slice_hits.size();

        double slice_total_charge = 0.0;
        double slice_neutrino_charge = 0.0;
        double slice_pattern_charge = 0.0;
        int slice_neutrino_hits = 0;
        int slice_pattern_hits = 0;

        for (const auto& hit : slice_hits) {
            slice_total_charge += hit->Integral();
            const auto& matched_mc_part = _assoc_mc_part->at(hit.key());
            const auto& matched_data = _assoc_mc_part->data(hit.key());
            double neutrino_charge = 0.0;
            double pattern_charge = 0.0;
            bool is_neutrino_hit = false;
            bool is_pattern_hit = false;
            for (size_t j = 0; j < matched_mc_part.size(); ++j) {
                int track_id = matched_mc_part[j]->TrackId();
                if (_neutrino_track_ids.count(track_id)) {
                    neutrino_charge += hit->Integral() * matched_data[j]->ideFraction;
                    if (matched_data[j]->isMaxIDE == 1) {
                        is_neutrino_hit = true;
                    }
                }
                if (_pattern_track_ids.count(track_id)) {
                    pattern_charge += hit->Integral() * matched_data[j]->ideFraction;
                    if (matched_data[j]->isMaxIDE == 1) {
                        is_pattern_hit = true;
                    }
                }
            }
            slice_neutrino_charge += neutrino_charge;
            slice_pattern_charge += pattern_charge;
            if (is_neutrino_hit) {
                ++slice_neutrino_hits;
            }
            if (is_pattern_hit) {
                ++slice_pattern_hits;
            }
        }

        _event_slice.total_charge = slice_total_charge;

        _event_slice.hit_completeness = (_total_neutrino_hits > 0) ? static_cast<float>(slice_neutrino_hits) / _total_neutrino_hits : 0.0f;
        _event_slice.hit_purity = (_event_slice.num_hits > 0) ? static_cast<float>(slice_neutrino_hits) / _event_slice.num_hits : 0.0f;
        _event_slice.charge_completeness = (_total_neutrino_charge > 0) ? slice_neutrino_charge / _total_neutrino_charge : 0.0f;
        _event_slice.charge_purity = (slice_total_charge > 0) ? slice_neutrino_charge / slice_total_charge : 0.0f;
        _event_slice.pattern_hit_purity = (_event_slice.num_hits > 0) ? static_cast<float>(slice_pattern_hits) / _event_slice.num_hits : 0.0f;
        _event_slice.pattern_hit_completeness = (_total_pattern_hits > 0) ? static_cast<float>(slice_pattern_hits) / _total_pattern_hits : 0.0f;
        _event_slice.pattern_charge_purity = (slice_total_charge > 0) ? slice_pattern_charge / slice_total_charge : 0.0f;
        _event_slice.pattern_charge_completeness = (_total_pattern_charge > 0) ? slice_pattern_charge / _total_pattern_charge : 0.0f;

        std::vector<art::Ptr<recob::SpacePoint>> slice_space_points;
        for (const auto& pfp : slice_pfp_v) {
            auto sp_v = pfp.get<recob::SpacePoint>();
            slice_space_points.insert(slice_space_points.end(), sp_v.begin(), sp_v.end());
        }

        std::vector<TVector3> positions;
        for (const auto& sp : slice_space_points) {
            const double* xyz = sp->XYZ();
            positions.emplace_back(xyz[0], xyz[1], xyz[2]);
        }

        if (positions.size() < 2) {
            _event_slice.pca_eigenvalues[0] = -1;
            _event_slice.pca_eigenvalues[1] = -1;
            _event_slice.pca_eigenvalues[2] = -1;
            _event_slice.pca_principal_dir[0] = 0;
            _event_slice.pca_principal_dir[1] = 0;
            _event_slice.pca_principal_dir[2] = 0;
        } else {
            TVector3 centroid(0, 0, 0);
            for (const auto& pos : positions) {
                centroid += pos;
            }
            centroid *= (1.0 / static_cast<double>(positions.size()));

            Eigen::MatrixXd centered(3, positions.size());
            for (size_t i = 0; i < positions.size(); ++i) {
                centered(0, i) = positions[i].X() - centroid.X();
                centered(1, i) = positions[i].Y() - centroid.Y();
                centered(2, i) = positions[i].Z() - centroid.Z();
            }

            Eigen::Matrix3d cov = (centered * centered.transpose()) / positions.size();

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(cov);
            if (eigensolver.info() == Eigen::Success) {
                Eigen::Vector3d eigenvalues = eigensolver.eigenvalues();
                Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();

                _event_slice.pca_eigenvalues[0] = eigenvalues(2);
                _event_slice.pca_eigenvalues[1] = eigenvalues(1);
                _event_slice.pca_eigenvalues[2] = eigenvalues(0);

                _event_slice.pca_principal_dir[0] = eigenvectors(0, 2);
                _event_slice.pca_principal_dir[1] = eigenvectors(1, 2);
                _event_slice.pca_principal_dir[2] = eigenvectors(2, 2);
            } else {
                _event_slice.pca_eigenvalues[0] = -1;
                _event_slice.pca_eigenvalues[1] = -1;
                _event_slice.pca_eigenvalues[2] = -1;
                _event_slice.pca_principal_dir[0] = 0;
                _event_slice.pca_principal_dir[1] = 0;
                _event_slice.pca_principal_dir[2] = 0;
            }
        }
    }

    void SliceAnalysis::setBranches(TTree* tree) {
        tree->Branch("slice_hit_completeness", &_event_slice.hit_completeness, "slice_hit_completeness/F");
        tree->Branch("slice_hit_purity", &_event_slice.hit_purity, "slice_hit_purity/F");
        tree->Branch("slice_charge_completeness", &_event_slice.charge_completeness, "slice_charge_completeness/F");
        tree->Branch("slice_charge_purity", &_event_slice.charge_purity, "slice_charge_purity/F");
        tree->Branch("pattern_hit_purity", &_event_slice.pattern_hit_purity, "pattern_hit_purity/F");
        tree->Branch("pattern_hit_completeness", &_event_slice.pattern_hit_completeness, "pattern_hit_completeness/F");
        tree->Branch("pattern_charge_purity", &_event_slice.pattern_charge_purity, "pattern_charge_purity/F");
        tree->Branch("pattern_charge_completeness", &_event_slice.pattern_charge_completeness, "pattern_charge_completeness/F");
        tree->Branch("slice_total_charge", &_event_slice.total_charge, "slice_total_charge/F");
        tree->Branch("pca_eigenvalues", _event_slice.pca_eigenvalues, "pca_eigenvalues[3]/F");
        tree->Branch("pca_principal_dir", _event_slice.pca_principal_dir, "pca_principal_dir[3]/F");
    }

    void SliceAnalysis::resetTTree(TTree* tree) {
        _event_slice = Slice();
    }

    DEFINE_ART_CLASS_TOOL(SliceAnalysis)
}

#endif