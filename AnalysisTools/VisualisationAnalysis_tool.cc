#ifndef ANALYSIS_VISUALISATION_CXX
#define ANALYSIS_VISUALISATION_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "TVector3.h"

#include "../CommonFuncs/BacktrackingFuncs.h"
#include "../CommonFuncs/TrackShowerScoreFuncs.h"
#include "../CommonFuncs/SpaceChargeCorrections.h"
#include "../CommonFuncs/Scatters.h"
#include "../CommonFuncs/PandoraFuncs.h"
#include "../CommonFuncs/TypeDefs.h"

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

    void analyzeEvent(art::Event const &e, bool fData) override;

    void analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

    void SaveTruth(art::Event const &e);

    void setBranches(TTree *_tree) override;

    void resetTTree(TTree *_tree) override;

private:
    std::vector<double> _all_hits_u_wire;
    std::vector<double> _all_hits_u_drift;
    std::vector<double> _all_hits_u_owner;
    std::vector<double> _all_hits_v_wire;
    std::vector<double> _all_hits_v_drift;
    std::vector<double> _all_hits_v_owner;
    std::vector<double> _all_hits_w_wire;
    std::vector<double> _all_hits_w_drift;
    std::vector<double> _all_hits_w_owner;

    std::vector<double> _slice_hits_u_wire;
    std::vector<double> _slice_hits_u_drift;
    std::vector<double> _slice_hits_v_wire;
    std::vector<double> _slice_hits_v_drift;
    std::vector<double> _slice_hits_w_wire;
    std::vector<double> _slice_hits_w_drift;

    art::Handle<std::vector<simb::MCParticle>> _mc_particle_handle; 
    std::vector<art::Ptr<simb::MCParticle>> _mc_particle_vector;
    lar_pandora::MCParticleMap _mc_particle_map;

    art::Handle<std::vector<recob::PFParticle>> _pf_particle_handle;
    std::vector<art::Ptr<recob::PFParticle>> _pf_particle_vector;
    lar_pandora::PFParticleMap _pf_particle_map;

    art::Handle<std::vector<recob::Hit>> _hit_handle;
    std::vector<art::Ptr<recob::Hit>> _hit_vector;

    art::Handle<std::vector<recob::Slice>> _slice_handle;
    std::vector<art::Ptr<recob::Slice>> _slice_vector;

    art::Handle<std::vector<recob::PFParticle>> _flash_match_pfp_handle;
    std::vector<art::Ptr<recob::PFParticle>> _flash_match_pfp_vector;
    std::vector<art::Ptr<recob::PFParticle>> _flash_nu_pfp_vector;

    std::map<int, int> _hits_to_track_id;
    std::map<int, std::vector<art::Ptr<recob::Hit>>> _track_id_to_hits;

    // Module labels for input tags
    std::string _MCParticleModuleLabel;
    std::string _PandoraModuleLabel;
    std::string _HitModuleLabel;
    std::string _BacktrackModuleLabel;
    std::string _FlashMatchModuleLabel;

    art::InputTag fCLSproducer;
};

VisualisationAnalysis::VisualisationAnalysis(const fhicl::ParameterSet &pset)
{
    _MCParticleModuleLabel = pset.get<std::string>("MCParticleModuleLabel");
    _PandoraModuleLabel = pset.get<std::string>("PandoraModuleLabel");
    _HitModuleLabel = pset.get<std::string>("HitModuleLabel");
    _BacktrackModuleLabel = pset.get<std::string>("BacktrackModuleLabel");
    _FlashMatchModuleLabel = pset.get<std::string>("FlashMatchModuleLabel");

    fCLSproducer = pset.get<art::InputTag>("CLSproducer");
}

void VisualisationAnalysis::configure(fhicl::ParameterSet const &pset)
{
}

void VisualisationAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
    if (!e.getByLabel(_MCParticleModuleLabel, _mc_particle_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any mc particles in event" << std::endl;
    art::fill_ptr_vector(_mc_particle_vector, _mc_particle_handle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(_mc_particle_vector, _mc_particle_map);

    if (!e.getByLabel(_PandoraModuleLabel, _pf_particle_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any pandora-slice pf particles in event" << std::endl;
    art::fill_ptr_vector(_pf_particle_vector, _pf_particle_handle);
    lar_pandora::LArPandoraHelper::BuildPFParticleMap(_pf_particle_vector, _pf_particle_map);

    if (!e.getByLabel(_HitModuleLabel, _hit_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any hits in event" << std::endl;
    art::fill_ptr_vector(_hit_vector, _hit_handle);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(_hit_handle, e, _BacktrackModuleLabel);

    if (!e.getByLabel(_PandoraModuleLabel, _slice_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any pandora slices in event" << std::endl;

    art::FindManyP<recob::Hit> hit_assoc = art::FindManyP<recob::Hit>(_slice_handle, e, _PandoraModuleLabel);
    art::fill_ptr_vector(_slice_vector, _slice_handle);

    art::FindManyP<recob::PFParticle> pf_part_slice_assoc = art::FindManyP<recob::PFParticle>(_slice_handle, e, _PandoraModuleLabel);
    //art::FindManyP<lardataobj::PFParticleMetadata> pf_part_metadata_assoc = art::FindManyP<lardataobj::PFParticleMetadata>(_pf_particle_handle, e, _PandoraModuleLabel);

    if (!e.getByLabel(_FlashMatchModuleLabel, _flash_match_pfp_handle))
        throw cet::exception("VisualisationAnalysis") << "failed to find any flash-matched pf particles" << std::endl;

    art::fill_ptr_vector(_flash_match_pfp_vector, _flash_match_pfp_handle);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(_flash_match_pfp_vector, _flash_nu_pfp_vector);

    for (unsigned int i_h = 0; i_h < _hit_vector.size(); i_h++)
    {
        const art::Ptr<recob::Hit> &hit = _hit_vector[i_h];
        const std::vector<art::Ptr<simb::MCParticle>> &matched_mc_part_vector = assoc_mc_part.at(hit.key());
        auto matched_data_vector = assoc_mc_part.data(hit.key());

        for (unsigned int i_p = 0; i_p < matched_mc_part_vector.size(); i_p++)
        {
            const art::Ptr<simb::MCParticle> &matched_mc_part = matched_mc_part_vector.at(i_p);
            auto matched_data = matched_data_vector.at(i_p);

            if (matched_data->isMaxIDE != 1)
                continue;

            const int track_id = matched_mc_part->TrackId();
            _hits_to_track_id[hit.key()] = track_id;
            _track_id_to_hits[track_id].push_back(hit);
        }
    }

    // Process hits in the event
    for (const art::Ptr<recob::Hit> &hit : _hit_vector)
    {
        common::PandoraView pandora_view = common::GetPandoraView(hit);
        TVector3 pandora_pos = common::GetPandoraHitPosition(e, hit, pandora_view);

        if (pandora_view == common::TPC_VIEW_U) {
            // Store hit information for the U-plane
            _all_hits_u_wire.push_back(pandora_pos.Z());
            _all_hits_u_drift.push_back(pandora_pos.X());
            _all_hits_u_owner.push_back((_hits_to_track_id.find(hit.key()) == _hits_to_track_id.end()) ?
                                        std::numeric_limits<double>::lowest() :
                                        _mc_particle_map.at(_hits_to_track_id.at(hit.key()))->PdgCode());
        }
        else if (pandora_view == common::TPC_VIEW_V) {
            // Store hit information for the V-plane
            _all_hits_v_wire.push_back(pandora_pos.Z());
            _all_hits_v_drift.push_back(pandora_pos.X());
            _all_hits_v_owner.push_back((_hits_to_track_id.find(hit.key()) == _hits_to_track_id.end()) ?
                                        std::numeric_limits<double>::lowest() :
                                        _mc_particle_map.at(_hits_to_track_id.at(hit.key()))->PdgCode());
        }
        else if (pandora_view == common::TPC_VIEW_W) {
            // Store hit information for the W-plane
            _all_hits_w_wire.push_back(pandora_pos.Z());
            _all_hits_w_drift.push_back(pandora_pos.X());
            _all_hits_w_owner.push_back((_hits_to_track_id.find(hit.key()) == _hits_to_track_id.end()) ?
                                        std::numeric_limits<double>::lowest() :
                                        _mc_particle_map.at(_hits_to_track_id.at(hit.key()))->PdgCode());
        }
    }
}

void VisualisationAnalysis::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                            proxy::withAssociated<recob::Hit>(fCLSproducer));

    for (const auto& pfp : slice_pfp_v)
    {
        if (pfp->IsPrimary())
            continue;

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
                _slice_hits_u_wire.push_back(pandora_pos.Z());
                _slice_hits_u_drift.push_back(pandora_pos.X());
            }
            else if (pandora_view == common::TPC_VIEW_V)
            {
                _slice_hits_v_wire.push_back(pandora_pos.Z());
                _slice_hits_v_drift.push_back(pandora_pos.X());
            }
            else if (pandora_view == common::TPC_VIEW_W)
            {
                _slice_hits_w_wire.push_back(pandora_pos.Z());
                _slice_hits_w_drift.push_back(pandora_pos.X());
            }
        }
    }
}


void VisualisationAnalysis::setBranches(TTree *_tree)
{
    _tree->Branch("all_hits_u_wire", &_all_hits_u_wire);
    _tree->Branch("all_hits_u_drift", &_all_hits_u_drift);
    _tree->Branch("all_hits_u_owner", &_all_hits_u_owner);
    _tree->Branch("all_hits_v_wire", &_all_hits_v_wire);
    _tree->Branch("all_hits_v_drift", &_all_hits_v_drift);
    _tree->Branch("all_hits_v_owner", &_all_hits_v_owner);
    _tree->Branch("all_hits_w_wire", &_all_hits_w_wire);
    _tree->Branch("all_hits_w_drift", &_all_hits_w_drift);
    _tree->Branch("all_hits_w_owner", &_all_hits_w_owner);

    _tree->Branch("slice_hits_u_wire", &_slice_hits_u_wire);
    _tree->Branch("slice_hits_u_drift", &_slice_hits_u_drift);
    _tree->Branch("slice_hits_v_wire", &_slice_hits_v_wire);
    _tree->Branch("slice_hits_v_drift", &_slice_hits_v_drift);
    _tree->Branch("slice_hits_w_wire", &_slice_hits_w_wire);
    _tree->Branch("slice_hits_w_drift", &_slice_hits_w_drift);
}

void VisualisationAnalysis::resetTTree(TTree *_tree)
{
    _all_hits_u_wire.clear();
    _all_hits_u_drift.clear();
    _all_hits_u_owner.clear();
    _all_hits_v_wire.clear();
    _all_hits_v_drift.clear();
    _all_hits_v_owner.clear();
    _all_hits_w_wire.clear();
    _all_hits_w_drift.clear();
    _all_hits_w_owner.clear();

    _slice_hits_u_wire.clear();
    _slice_hits_u_drift.clear();
    _slice_hits_v_wire.clear();
    _slice_hits_v_drift.clear();
    _slice_hits_w_wire.clear();
    _slice_hits_w_drift.clear();
}

DEFINE_ART_CLASS_TOOL(VisualisationAnalysis)
} 

#endif
