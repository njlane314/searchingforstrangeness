#ifndef PATTERN_ANALYSIS_CXX
#define PATTERN_ANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "TVector3.h"
#include "CommonFunctions/Backtracking.h"
#include "CommonFunctions/Scores.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Types.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"
#include "EventClassifier.h"
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <memory>

namespace analysis
{
    class PatternAnalysis : public AnalysisToolBase {
    public:
        PatternAnalysis(const fhicl::ParameterSet &pset);
        ~PatternAnalysis() = default;

        void configure(fhicl::ParameterSet const &pset) override;
        void analyseEvent(art::Event const &e, bool is_data) override;
        void analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected) override;
        void setBranches(TTree *_tree) override;
        void resetTTree(TTree *_tree) override;

    private:
        void buildHitsToTrackMap(const std::vector<art::Ptr<recob::Hit>>& hit_vector,
                                 const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& assoc_mc_part,
                                 const lar_pandora::MCParticleMap& mc_particle_map,
                                 std::map<int, int>& hits_to_track_map,
                                 std::map<int, std::vector<art::Ptr<recob::Hit>>>& track_to_hits_map);

        void storeTrueHits(const std::vector<art::Ptr<recob::Hit>>& hit_vector,
                           const std::map<int, int>& hits_to_track_map,
                           const lar_pandora::MCParticleMap& mc_particle_map,
                           art::Event const &e);

        std::vector<art::Ptr<recob::Hit>> getPfpHits(const common::ProxyPfpElem_t& pfp,
                                                     const common::ProxyClusColl_t& clus_proxy);

        void getViewHits(const std::vector<art::Ptr<recob::Hit>>& hit_v,
                         art::Event const &e,
                         std::vector<float>& wire_u, std::vector<float>& drift_u,
                         std::vector<float>& wire_v, std::vector<float>& drift_v,
                         std::vector<float>& wire_w, std::vector<float>& drift_w);

        std::unique_ptr<signature::EventClassifier> _classifier;
        std::vector<std::set<size_t>> _signature_true_hits;
        std::vector<int> _signature_types;
        std::vector<float> _signature_completeness;
        std::vector<float> _signature_purity;

        std::vector<float> _true_hits_u_wire;
        std::vector<float> _true_hits_u_drift;
        std::vector<float> _true_hits_u_owner;
        std::vector<float> _true_hits_v_wire;
        std::vector<float> _true_hits_v_drift;
        std::vector<float> _true_hits_v_owner;
        std::vector<float> _true_hits_w_wire;
        std::vector<float> _true_hits_w_drift;
        std::vector<float> _true_hits_w_owner;

        std::vector<std::vector<float>> _slice_hits_u_wire;
        std::vector<std::vector<float>> _slice_hits_u_drift;
        std::vector<std::vector<float>> _slice_hits_v_wire;
        std::vector<std::vector<float>> _slice_hits_v_drift;
        std::vector<std::vector<float>> _slice_hits_w_wire;
        std::vector<std::vector<float>> _slice_hits_w_drift;

        art::InputTag _MCTproducer;
        art::InputTag _PFPproducer;
        art::InputTag _Hproducer;
        art::InputTag _BacktrackTag;
        art::InputTag _FMproducer;
        art::InputTag _CLSproducer;
    };

    PatternAnalysis::PatternAnalysis(const fhicl::ParameterSet &pset)
        : _MCTproducer(pset.get<art::InputTag>("MCTproducer", "largeant")),
        _PFPproducer(pset.get<art::InputTag>("PFPproducer", "pandoraPatRec:allOutcomes")),
        _Hproducer(pset.get<art::InputTag>("Hproducer", "gaushit")),
        _BacktrackTag(pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")),
        _FMproducer(pset.get<art::InputTag>("FMproducer", "pandora")),
        _CLSproducer(pset.get<art::InputTag>("CLSproducer", "pandora")),
        _classifier(std::make_unique<signature::EventClassifier>(pset.get<fhicl::ParameterSet>("EventClassifier")))
    {}

    void PatternAnalysis::configure(fhicl::ParameterSet const &pset) {}

    void PatternAnalysis::buildHitsToTrackMap(
        const std::vector<art::Ptr<recob::Hit>>& hit_vector,
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& assoc_mc_part,
        const lar_pandora::MCParticleMap& mc_particle_map,
        std::map<int, int>& hits_to_track_map,
        std::map<int, std::vector<art::Ptr<recob::Hit>>>& track_to_hits_map) {
        for (const auto& hit : hit_vector) {
            const auto& matched_mc_part_vector = assoc_mc_part.at(hit.key());
            const auto matched_data_vector = assoc_mc_part.data(hit.key());
            for (size_t i_p = 0; i_p < matched_mc_part_vector.size(); ++i_p) {
                const auto& matched_mc_part = matched_mc_part_vector[i_p];
                const auto matched_data = matched_data_vector[i_p];
                if (matched_data->isMaxIDE != 1) continue;
                int track_idx = common::isParticleElectromagnetic(matched_mc_part)
                                    ? common::getLeadElectromagneticTrack(matched_mc_part, mc_particle_map)
                                    : matched_mc_part->TrackId();
                hits_to_track_map[hit.key()] = track_idx;
                track_to_hits_map[track_idx].push_back(hit);
            }
        }
    }

    void PatternAnalysis::storeTrueHits(
        const std::vector<art::Ptr<recob::Hit>>& hit_vector,
        const std::map<int, int>& hits_to_track_map,
        const lar_pandora::MCParticleMap& mc_particle_map,
        art::Event const &e) {
        for (const auto& hit : hit_vector) {
            common::PandoraView pandora_view = common::GetPandoraView(hit);
            TVector3 pandora_pos = common::GetPandoraHitPosition(e, hit, pandora_view);
            auto hit_to_track_it = hits_to_track_map.find(hit.key());
            if (hit_to_track_it == hits_to_track_map.end()) continue;
            int owner_pdg_code = mc_particle_map.at(hit_to_track_it->second)->PdgCode();

            if (pandora_view == common::TPC_VIEW_U) {
                _true_hits_u_wire.push_back(pandora_pos.Z());
                _true_hits_u_drift.push_back(pandora_pos.X());
                _true_hits_u_owner.push_back(owner_pdg_code);
            } else if (pandora_view == common::TPC_VIEW_V) {
                _true_hits_v_wire.push_back(pandora_pos.Z());
                _true_hits_v_drift.push_back(pandora_pos.X());
                _true_hits_v_owner.push_back(owner_pdg_code);
            } else if (pandora_view == common::TPC_VIEW_W) {
                _true_hits_w_wire.push_back(pandora_pos.Z());
                _true_hits_w_drift.push_back(pandora_pos.X());
                _true_hits_w_owner.push_back(owner_pdg_code);
            }
        }
    }

    void PatternAnalysis::analyseEvent(art::Event const &e, bool is_data) {
        if (is_data) return;

        art::Handle<std::vector<simb::MCParticle>> mc_particle_handle;
        std::vector<art::Ptr<simb::MCParticle>> mc_particle_vector;
        lar_pandora::MCParticleMap mc_particle_map;
        if (e.getByLabel(_MCTproducer, mc_particle_handle)) {
            art::fill_ptr_vector(mc_particle_vector, mc_particle_handle);
            lar_pandora::LArPandoraHelper::BuildMCParticleMap(mc_particle_vector, mc_particle_map);
        } else return;

        art::Handle<std::vector<recob::PFParticle>> pf_particle_handle;
        if (!e.getByLabel(_PFPproducer, pf_particle_handle)) return;

        art::Handle<std::vector<recob::Hit>> hit_handle;
        std::vector<art::Ptr<recob::Hit>> hit_vector;
        if (e.getByLabel(_Hproducer, hit_handle)) {
            art::fill_ptr_vector(hit_vector, hit_handle);
            art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part(hit_handle, e, _BacktrackTag);

            std::map<int, int> hits_to_track_map;
            std::map<int, std::vector<art::Ptr<recob::Hit>>> track_to_hits_map;
            this->buildHitsToTrackMap(hit_vector, assoc_mc_part, mc_particle_map, hits_to_track_map, track_to_hits_map);
            this->storeTrueHits(hit_vector, hits_to_track_map, mc_particle_map, e);

            const auto& pattern = _classifier->getPattern(e);
            _signature_true_hits.clear();
            for (const auto& entry : pattern) {
                const auto& signature = entry.second;
                std::set<int> signature_track_ids;
                for (const auto& mcp : signature) {
                    int track_idx = common::isParticleElectromagnetic(mcp)
                                        ? common::getLeadElectromagneticTrack(mcp, mc_particle_map)
                                        : mcp->TrackId();
                    if (track_to_hits_map.count(track_idx)) {
                        signature_track_ids.insert(track_idx);
                    }
                }
                std::set<size_t> true_hits;
                for (int track_id : signature_track_ids) {
                    for (const auto& hit : track_to_hits_map[track_id]) {
                        true_hits.insert(hit.key());
                    }
                }
                _signature_true_hits.push_back(true_hits);
            }
        }
    }

    std::vector<art::Ptr<recob::Hit>> PatternAnalysis::getPfpHits(
        const common::ProxyPfpElem_t& pfp,
        const common::ProxyClusColl_t& clus_proxy) {
        std::vector<art::Ptr<recob::Hit>> hit_v;
        auto clus_pxy_v = pfp.get<recob::Cluster>();
        for (const auto& ass_clus : clus_pxy_v) {
            const auto& clus = clus_proxy[ass_clus.key()];
            auto clus_hit_v = clus.get<recob::Hit>();
            hit_v.insert(hit_v.end(), clus_hit_v.begin(), clus_hit_v.end());
        }
        return hit_v;
    }

    void PatternAnalysis::getViewHits(
        const std::vector<art::Ptr<recob::Hit>>& hit_v,
        art::Event const &e,
        std::vector<float>& wire_u, std::vector<float>& drift_u,
        std::vector<float>& wire_v, std::vector<float>& drift_v,
        std::vector<float>& wire_w, std::vector<float>& drift_w) {
        for (const auto& hit : hit_v) {
            common::PandoraView pandora_view = common::GetPandoraView(hit);
            TVector3 pandora_pos = common::GetPandoraHitPosition(e, hit, pandora_view);
            if (pandora_view == common::TPC_VIEW_U) {
                wire_u.push_back(pandora_pos.Z());
                drift_u.push_back(pandora_pos.X());
            } else if (pandora_view == common::TPC_VIEW_V) {
                wire_v.push_back(pandora_pos.Z());
                drift_v.push_back(pandora_pos.X());
            } else if (pandora_view == common::TPC_VIEW_W) {
                wire_w.push_back(pandora_pos.Z());
                drift_w.push_back(pandora_pos.X());
            }
        }
    }

    void PatternAnalysis::analyseSlice(
        art::Event const &e,
        std::vector<common::ProxyPfpElem_t> &slice_pfp_v,
        bool is_data,
        bool selected) {
        if (is_data) return;

        common::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(
            e, _CLSproducer, proxy::withAssociated<recob::Hit>(_CLSproducer));

        for (const auto& pfp : slice_pfp_v) {
            if (pfp->IsPrimary()) continue;

            std::vector<float> pfp_wire_u, pfp_drift_u, pfp_wire_v, pfp_drift_v, pfp_wire_w, pfp_drift_w;
            auto hit_v = this->getPfpHits(pfp, clus_proxy);
            this->getViewHits(hit_v, e, pfp_wire_u, pfp_drift_u, pfp_wire_v, pfp_drift_v, pfp_wire_w, pfp_drift_w);

            _slice_hits_u_wire.push_back(pfp_wire_u);
            _slice_hits_u_drift.push_back(pfp_drift_u);
            _slice_hits_v_wire.push_back(pfp_wire_v);
            _slice_hits_v_drift.push_back(pfp_drift_v);
            _slice_hits_w_wire.push_back(pfp_wire_w);
            _slice_hits_w_drift.push_back(pfp_drift_w);
        }

        if (selected) {
            std::vector<std::set<size_t>> pfp_hits_v;
            for (const auto& pfp : slice_pfp_v) {
                if (pfp->IsPrimary()) continue;
                auto hit_v = this->getPfpHits(pfp, clus_proxy);
                std::set<size_t> pfp_hits;
                for (const auto& hit : hit_v) {
                    pfp_hits.insert(hit.key());
                }
                pfp_hits_v.push_back(pfp_hits);
            }

            const auto& pattern = _classifier->getPattern(e);
            for (size_t s = 0; s < pattern.size(); ++s) {
                const auto& true_hits = _signature_true_hits[s];
                if (true_hits.empty()) continue;

                int max_shared = 0;
                size_t best_pfp_idx = 0;
                for (size_t i = 0; i < pfp_hits_v.size(); ++i) {
                    const auto& pfp_hits = pfp_hits_v[i];
                    std::vector<size_t> intersection;
                    std::set_intersection(true_hits.begin(), true_hits.end(), pfp_hits.begin(), pfp_hits.end(), std::back_inserter(intersection));
                    int shared = intersection.size();
                    if (shared > max_shared) {
                        max_shared = shared;
                        best_pfp_idx = i;
                    }
                }

                float completeness = (true_hits.size() > 0) ? static_cast<float>(max_shared) / true_hits.size() : 0.0f;
                float purity = (pfp_hits_v.size() > 0 && pfp_hits_v[best_pfp_idx].size() > 0) ? static_cast<float>(max_shared) / pfp_hits_v[best_pfp_idx].size() : 0.0f;
                _signature_types.push_back(static_cast<int>(pattern[s].first));
                _signature_completeness.push_back(completeness);
                _signature_purity.push_back(purity);
            }
        }
    }

    void PatternAnalysis::setBranches(TTree *_tree) {
        _tree->Branch("true_hits_u_wire", &_true_hits_u_wire);
        _tree->Branch("true_hits_u_drift", &_true_hits_u_drift);
        _tree->Branch("true_hits_u_owner", &_true_hits_u_owner);
        _tree->Branch("true_hits_v_wire", &_true_hits_v_wire);
        _tree->Branch("true_hits_v_drift", &_true_hits_v_drift);
        _tree->Branch("true_hits_v_owner", &_true_hits_v_owner);
        _tree->Branch("true_hits_w_wire", &_true_hits_w_wire);
        _tree->Branch("true_hits_w_drift", &_true_hits_w_drift);
        _tree->Branch("true_hits_w_owner", &_true_hits_w_owner);

        _tree->Branch("slice_hits_u_wire", &_slice_hits_u_wire);
        _tree->Branch("slice_hits_u_drift", &_slice_hits_u_drift);
        _tree->Branch("slice_hits_v_wire", &_slice_hits_v_wire);
        _tree->Branch("slice_hits_v_drift", &_slice_hits_v_drift);
        _tree->Branch("slice_hits_w_wire", &_slice_hits_w_wire);
        _tree->Branch("slice_hits_w_drift", &_slice_hits_w_drift);

        _tree->Branch("signature_types", &_signature_types);
        _tree->Branch("signature_completeness", &_signature_completeness);
        _tree->Branch("signature_purity", &_signature_purity);
    }

    void PatternAnalysis::resetTTree(TTree *_tree) {
        _true_hits_u_wire.clear();
        _true_hits_u_drift.clear();
        _true_hits_u_owner.clear();
        _true_hits_v_wire.clear();
        _true_hits_v_drift.clear();
        _true_hits_v_owner.clear();
        _true_hits_w_wire.clear();
        _true_hits_w_drift.clear();
        _true_hits_w_owner.clear();

        _slice_hits_u_wire.clear();
        _slice_hits_u_drift.clear();
        _slice_hits_v_wire.clear();
        _slice_hits_v_drift.clear();
        _slice_hits_w_wire.clear();
        _slice_hits_w_drift.clear();

        _signature_true_hits.clear();
        _signature_types.clear();
        _signature_completeness.clear();
        _signature_purity.clear();
    }

    DEFINE_ART_CLASS_TOOL(PatternAnalysis)
}

#endif