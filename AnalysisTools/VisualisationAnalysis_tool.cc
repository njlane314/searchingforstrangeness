#ifndef ANALYSIS_VISUALISATION_CXX
#define ANALYSIS_VISUALISATION_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "TVector3.h"

#include "../CommonFunctions/BacktrackingFuncs.h"
#include "../CommonFunctions/TrackShowerScoreFuncs.h"
#include "../CommonFunctions/SpaceChargeCorrections.h"
#include "../CommonFunctions/Scatters.h"
#include "../CommonFunctions/PandoraFuncs.h"
#include "../CommonFunctions/TypeDefs.h"

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

class VisualisationAnalysis : public AnalysisToolBase
{

public:
  
    VisualisationAnalysis(const fhicl::ParameterSet &pset);

    ~VisualisationAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);

    void analyzeEvent(art::Event const &e, bool is_data) override;

    void analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected) override;

    void SaveTruth(art::Event const &e);

    void setBranches(TTree *_tree) override;

    void resetTTree(TTree *_tree) override;

private:

    bool isParticleElectromagnetic(const art::Ptr<simb::MCParticle> &mc_part);
    int getLeadElectromagneticTrack(const art::Ptr<simb::MCParticle> &mc_part, const lar_pandora::MCParticleMap &mc_particle_map);

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

    std::string _SimulationModuleLabel;
    std::string _PandoraModuleLabel;
    std::string _HitModuleLabel;
    std::string _BacktrackModuleLabel;
    std::string _FlashMatchModuleLabel;

    art::InputTag _CLSproducer;
};

VisualisationAnalysis::VisualisationAnalysis(const fhicl::ParameterSet &pset)
{
    _SimulationModuleLabel = pset.get<std::string>("MCParticleModuleLabel");
    _PandoraModuleLabel = pset.get<std::string>("PandoraModuleLabel");
    _HitModuleLabel = pset.get<std::string>("HitModuleLabel");
    _BacktrackModuleLabel = pset.get<std::string>("BacktrackModuleLabel");
    _FlashMatchModuleLabel = pset.get<std::string>("FlashMatchModuleLabel");

    _CLSproducer = pset.get<art::InputTag>("CLSproducer");
}

void VisualisationAnalysis::configure(fhicl::ParameterSet const &pset)
{
}

void VisualisationAnalysis::analyzeEvent(art::Event const &e, bool is_data)
{
     art::Handle<std::vector<simb::MCParticle>> mc_particle_handle; 
    std::vector<art::Ptr<simb::MCParticle>> mc_particle_vector;
    lar_pandora::MCParticleMap mc_particle_map;

    if (!e.getByLabel(_SimulationModuleLabel, mc_particle_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any mc particles in event" << std::endl;
    art::fill_ptr_vector(mc_particle_vector, mc_particle_handle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mc_particle_vector, mc_particle_map);

    art::Handle<std::vector<recob::PFParticle>> pf_particle_handle;
    std::vector<art::Ptr<recob::PFParticle>> pf_particle_vector;
    lar_pandora::PFParticleMap pf_particle_map;
    
    if (!e.getByLabel(_PandoraModuleLabel, pf_particle_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any pandora-slice pf particles in event" << std::endl;
    art::fill_ptr_vector(pf_particle_vector, pf_particle_handle);
    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pf_particle_vector, pf_particle_map);

    art::Handle<std::vector<recob::Hit>> hit_handle;
    std::vector<art::Ptr<recob::Hit>> hit_vector;
    
    if (!e.getByLabel(_HitModuleLabel, hit_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any hits in event" << std::endl;
    art::fill_ptr_vector(hit_vector, hit_handle);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hit_handle, e, _BacktrackModuleLabel);

    art::Handle<std::vector<recob::Slice>> slice_handle;
    std::vector<art::Ptr<recob::Slice>> slice_vector;

    if (!e.getByLabel(_PandoraModuleLabel, slice_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any pandora slices in event" << std::endl;
    art::FindManyP<recob::Hit> hit_assoc = art::FindManyP<recob::Hit>(slice_handle, e, _PandoraModuleLabel);
    art::fill_ptr_vector(slice_vector, slice_handle);
    art::FindManyP<recob::PFParticle> pf_part_slice_assoc = art::FindManyP<recob::PFParticle>(slice_handle, e, _PandoraModuleLabel);
    //art::FindManyP<lardataobj::PFParticleMetadata> pf_part_metadata_assoc = art::FindManyP<lardataobj::PFParticleMetadata>(pf_particle_handle, e, _PandoraModuleLabel);

    art::Handle<std::vector<recob::PFParticle>> flash_match_pf_particle_handle;
    std::vector<art::Ptr<recob::PFParticle>> flash_match_pf_particle_vector;
    std::vector<art::Ptr<recob::PFParticle>> flash_nu_pf_particle_vector;

    if (!e.getByLabel(_FlashMatchModuleLabel, flash_match_pf_particle_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any flash-matched pf particles" << std::endl;
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

            const int track_idx = isParticleElectromagnetic(matched_mc_part) ? getLeadElectromagneticTrack(matched_mc_part, mc_particle_map) : matched_mc_part->TrackId();

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
}

void VisualisationAnalysis::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected)
{
    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, _CLSproducer,
                                                                                            proxy::withAssociated<recob::Hit>(_CLSproducer));

    std::cout << "Size of slice " << slice_pfp_v.size() << std::endl;

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

void VisualisationAnalysis::setBranches(TTree *_tree)
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

void VisualisationAnalysis::resetTTree(TTree *_tree)
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

bool VisualisationAnalysis::isParticleElectromagnetic(const art::Ptr<simb::MCParticle> &mc_part)
{
    return ((std::abs(mc_part->PdgCode() == 11) || (mc_part->PdgCode() == 22)));
}

int VisualisationAnalysis::getLeadElectromagneticTrack(const art::Ptr<simb::MCParticle> &mc_part, const lar_pandora::MCParticleMap &mc_particle_map)
{
    int track_idx = mc_part->TrackId();
    art::Ptr<simb::MCParticle> mother_mc_part = mc_part;

    do 
    {
        track_idx = mother_mc_part->TrackId();
        const int mother_idx = mother_mc_part->Mother();

        if (mc_particle_map.find(mother_idx) == mc_particle_map.end())
            break;

        mother_mc_part = mc_particle_map.at(mother_idx);
    } 
    while (isParticleElectromagnetic(mother_mc_part));

    return track_idx;
}


DEFINE_ART_CLASS_TOOL(VisualisationAnalysis)
} 

#endif