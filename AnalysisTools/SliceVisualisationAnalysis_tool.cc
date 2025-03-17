#ifndef ANALYSIS_SLICEVISUALISATION_CXX
#define ANALYSIS_SLICEVISUALISATION_CXX

#include <iostream>
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

namespace analysis
{

class SliceVisualisationAnalysis : public AnalysisToolBase
{

public:
    SliceVisualisationAnalysis(const fhicl::ParameterSet &pset);
    ~SliceVisualisationAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);
    void analyseEvent(art::Event const &e, bool is_data) override;
    void analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected) override;
    void SaveTruth(art::Event const &e);
    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

private:
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

SliceVisualisationAnalysis::SliceVisualisationAnalysis(const fhicl::ParameterSet &pset)
{
    _MCTproducer = pset.get<art::InputTag>("MCTproducer", "largeant");
    _PFPproducer = pset.get<art::InputTag>("PFPproducer", "pandoraPatRec:allOutcomes");
    _Hproducer = pset.get<art::InputTag>("Hproducer", "gaushit");
    _BacktrackTag = pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
    _FMproducer = pset.get<art::InputTag>("FMproducer", "pandora");
    _CLSproducer = pset.get<art::InputTag>("CLSproducer", "pandora");
}

void SliceVisualisationAnalysis::configure(fhicl::ParameterSet const &pset)
{
}

void SliceVisualisationAnalysis::analyseEvent(art::Event const &e, bool is_data)
{
    std::cout << "Analysing event in SliceVisualisation..." << std::endl;

    art::Handle<std::vector<simb::MCParticle>> mc_particle_handle; 
    std::vector<art::Ptr<simb::MCParticle>> mc_particle_vector;
    lar_pandora::MCParticleMap mc_particle_map;

    if (!e.getByLabel(_MCTproducer, mc_particle_handle))
        throw cet::exception("SliceVisualisationAnalysis") << "failed to find any mc particles in event" << std::endl;
    art::fill_ptr_vector(mc_particle_vector, mc_particle_handle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mc_particle_vector, mc_particle_map);

    art::Handle<std::vector<recob::PFParticle>> pf_particle_handle;
    std::vector<art::Ptr<recob::PFParticle>> pf_particle_vector;
    lar_pandora::PFParticleMap pf_particle_map;
    
    if (!e.getByLabel(_PFPproducer, pf_particle_handle))
        throw cet::exception("SliceVisualisationAnalysis") << "failed to find any pandora-slice pf particles in event" << std::endl;
    art::fill_ptr_vector(pf_particle_vector, pf_particle_handle);
    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pf_particle_vector, pf_particle_map);

    art::Handle<std::vector<recob::Hit>> hit_handle;
    std::vector<art::Ptr<recob::Hit>> hit_vector;
    
    if (!e.getByLabel(_Hproducer, hit_handle))
        throw cet::exception("SliceVisualisationAnalysis") << "failed to find any hits in event" << std::endl;
    art::fill_ptr_vector(hit_vector, hit_handle);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hit_handle, e, _BacktrackTag);

    art::Handle<std::vector<recob::Slice>> slice_handle;
    std::vector<art::Ptr<recob::Slice>> slice_vector;

    if (!e.getByLabel(_PFPproducer, slice_handle))
        throw cet::exception("SliceVisualisationAnalysis") << "failed to find any pandora slices in event" << std::endl;
    art::FindManyP<recob::Hit> hit_assoc = art::FindManyP<recob::Hit>(slice_handle, e, _PFPproducer);
    art::fill_ptr_vector(slice_vector, slice_handle);
    art::FindManyP<recob::PFParticle> pf_part_slice_assoc = art::FindManyP<recob::PFParticle>(slice_handle, e, _PFPproducer);

    art::Handle<std::vector<recob::PFParticle>> flash_match_pf_particle_handle;
    std::vector<art::Ptr<recob::PFParticle>> flash_match_pf_particle_vector;
    std::vector<art::Ptr<recob::PFParticle>> flash_nu_pf_particle_vector;

    if (!e.getByLabel(_FMproducer, flash_match_pf_particle_handle))
        throw cet::exception("SliceVisualisationAnalysis") << "failed to find any flash-matched pf particles" << std::endl;
    art::fill_ptr_vector(flash_match_pf_particle_vector, flash_match_pf_particle_handle);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(flash_match_pf_particle_vector, flash_nu_pf_particle_vector);

    std::map<int, int> hits_to_track_map;
    std::map<int, std::vector<art::Ptr<recob::Hit>>> track_to_hits_map;

    for (unsigned int i_h = 0; i_h < hit_vector.size(); i_h++)
    {
        const art::Ptr<recob::Hit> &hit = hit_vector[i_h];
        const std::vector<art::Ptr<simb::MCParticle>> &matched_mc_part_vector = assoc_mc_part.at(hit.key());
        auto matched_data_vector = assoc_mc_part.data(hit.key());

        for (unsigned int i_p = 0; i_p < matched_mc_part_vector.size(); i_p++)
        {
            const art::Ptr<simb::MCParticle> &matched_mc_part = matched_mc_part_vector.at(i_p);
            auto matched_data = matched_data_vector.at(i_p);

            if (matched_data->isMaxIDE != 1)
                continue;

            const int track_idx = common::isParticleElectromagnetic(matched_mc_part) ? common::getLeadElectromagneticTrack(matched_mc_part, mc_particle_map) : matched_mc_part->TrackId();

            hits_to_track_map[hit.key()] = track_idx;
            track_to_hits_map[track_idx].push_back(hit);
        }
    }

    // Process hits in the event
    for (const art::Ptr<recob::Hit> &hit : hit_vector)
    {
        common::PandoraView pandora_view = common::GetPandoraView(hit);
        TVector3 pandora_pos = common::GetPandoraHitPosition(e, hit, pandora_view);

        auto hit_to_track_it = hits_to_track_map.find(hit.key());
        if (hit_to_track_it == hits_to_track_map.end()) {
            continue; 
        }

        int owner_pdg_code = mc_particle_map.at(hit_to_track_it->second)->PdgCode();

        if (pandora_view == common::TPC_VIEW_U) {
            // Store hit information for the U-plane
            _true_hits_u_wire.push_back(pandora_pos.Z());
            _true_hits_u_drift.push_back(pandora_pos.X());
            _true_hits_u_owner.push_back(owner_pdg_code);
        }
        else if (pandora_view == common::TPC_VIEW_V) {
            // Store hit information for the V-plane
            _true_hits_v_wire.push_back(pandora_pos.Z());
            _true_hits_v_drift.push_back(pandora_pos.X());
            _true_hits_v_owner.push_back(owner_pdg_code);
        }
        else if (pandora_view == common::TPC_VIEW_W) {
            // Store hit information for the W-plane
            _true_hits_w_wire.push_back(pandora_pos.Z());
            _true_hits_w_drift.push_back(pandora_pos.X());
            _true_hits_w_owner.push_back(owner_pdg_code);
        }
    }

    std::cout << "Finished analysing event in SliceVisualisation!" << std::endl;
}

void SliceVisualisationAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected)
{
    std::cout << "Analysisng slice in SliceVisualisation..." << std::endl;
    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, _CLSproducer,
                                                                                            proxy::withAssociated<recob::Hit>(_CLSproducer));

    for (const auto& pfp : slice_pfp_v)
    {
        if (pfp->IsPrimary())
            continue;

        std::vector<float> pfp_slice_hits_u_wire;
        std::vector<float> pfp_slice_hits_u_drift;
        std::vector<float> pfp_slice_hits_v_wire;
        std::vector<float> pfp_slice_hits_v_drift;
        std::vector<float> pfp_slice_hits_w_wire;
        std::vector<float> pfp_slice_hits_w_drift;

        std::vector<art::Ptr<recob::Hit>> hit_v;
        auto clus_pxy_v = pfp.get<recob::Cluster>();

        for (auto ass_clus : clus_pxy_v)
        {
            const auto &clus = clus_proxy[ass_clus.key()];
            auto clus_hit_v = clus.get<recob::Hit>();
            
            for (size_t h = 0; h < clus_hit_v.size(); h++)
            {
                hit_v.push_back(clus_hit_v[h]);
            }
        }

        for (const auto& hit : hit_v)
        {
            common::PandoraView pandora_view = common::GetPandoraView(hit);
            TVector3 pandora_pos = common::GetPandoraHitPosition(e, hit, pandora_view);

            if (pandora_view == common::TPC_VIEW_U)
            {
                pfp_slice_hits_u_wire.push_back(pandora_pos.Z());
                pfp_slice_hits_u_drift.push_back(pandora_pos.X());
            }
            else if (pandora_view == common::TPC_VIEW_V)
            {
                pfp_slice_hits_v_wire.push_back(pandora_pos.Z());
                pfp_slice_hits_v_drift.push_back(pandora_pos.X());
            }
            else if (pandora_view == common::TPC_VIEW_W)
            {
                pfp_slice_hits_w_wire.push_back(pandora_pos.Z());
                pfp_slice_hits_w_drift.push_back(pandora_pos.X());
            }
        }

        _slice_hits_u_wire.push_back(pfp_slice_hits_u_wire);
        _slice_hits_u_drift.push_back(pfp_slice_hits_u_drift);
        _slice_hits_v_wire.push_back(pfp_slice_hits_v_wire);
        _slice_hits_v_drift.push_back(pfp_slice_hits_v_drift);
        _slice_hits_w_wire.push_back(pfp_slice_hits_w_wire);
        _slice_hits_w_drift.push_back(pfp_slice_hits_w_drift);
    }    
}

void SliceVisualisationAnalysis::setBranches(TTree *_tree)
{
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
}

void SliceVisualisationAnalysis::resetTTree(TTree *_tree)
{
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
}

DEFINE_ART_CLASS_TOOL(SliceVisualisationAnalysis)
} 

#endif