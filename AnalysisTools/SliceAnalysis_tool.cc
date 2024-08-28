#ifndef ANALYSIS_SLICE_CXX
#define ANALYSIS_SLICE_CXX

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

class SliceAnalysis : public AnalysisToolBase
{

public:
  
    SliceAnalysis(const fhicl::ParameterSet &pset);

    ~SliceAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);

    void analyzeEvent(art::Event const &e, bool fData) override;

    void analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

    void SaveTruth(art::Event const &e);

    void setBranches(TTree *_tree) override;

    void resetTTree(TTree *_tree) override;

private:

    bool isParticleElectromagnetic(const art::Ptr<simb::MCParticle> &mc_part);
    int getLeadElectromagneticTrack(const art::Ptr<simb::MCParticle> &mc_part);

    std::string _MCParticleModuleLabel;
    std::string _PandoraModuleLabel;
    std::string _HitModuleLabel;
    std::string _BacktrackModuleLabel;
    std::string _FlashMatchModuleLabel;

    art::Handle<std::vector<simb::MCParticle>> _mc_particle_handle;
    std::vector<art::Ptr<simb::MCParticle>> _mc_particle_vector;

    art::Handle<std::vector<recob::PFParticle>> _pf_particle_handle;
    std::vector<art::Ptr<recob::PFParticle>> _pf_particle_vector;

    art::Handle<std::vector<recob::Hit>> _hit_handle;
    std::vector<art::Ptr<recob::Hit>> _hit_vector;

    art::Handle<std::vector<recob::Slice>> _slice_handle; 
    std::vector<art::Ptr<recob::Slice>> _slice_vector;

    art::Handle<std::vector<recob::PFParticle>> _flash_match_pfp_handle;
    std::vector<art::Ptr<recob::PFParticle>> _flash_match_pfp_vector;
    std::vector<art::Ptr<recob::PFParticle>> _flash_nu_pfp_vector;

    lar_pandora::PFParticleMap _pf_particle_map;
    lar_pandora::MCParticleMap _mc_particle_map; 
    std::map<int, int> _hit_to_track_id;
    std::map<int, std::vector<art::Ptr<recob::Hit>>> _track_id_to_hits;

    int _true_nu_slice_id;
    int _pandora_nu_slice_id;
    int _flash_match_nu_slice_id;

    float _true_slice_completeness;
    float _true_slice_purity;
    float _pandora_slice_completeness;
    float _pandora_slice_purity;
    float _flash_slice_completeness;
    float _flash_slice_purity;

    float _true_reco_nu_vtx_x;
    float _true_reco_nu_vtx_y;
    float _true_reco_nu_vtx_z;
    float _pandora_reco_nu_vtx_x;
    float _pandora_reco_nu_vtx_y;
    float _pandora_reco_nu_vtx_z;
    float _flash_reco_nu_vtx_x;
    float _flash_reco_nu_vtx_y;
    float _flash_reco_nu_vtx_z;
};

SliceAnalysis::SliceAnalysis(const fhicl::ParameterSet &pset)
{
    _MCParticleModuleLabel = pset.get<std::string>("MCParticleModuleLabel", "largeant");
    _PandoraModuleLabel = pset.get<std::string>("PandoraModuleLabel", "pandora");
    _HitModuleLabel = pset.get<std::string>("HitModuleLabel", "gaushit");
    _BacktrackModuleLabel = pset.get<std::string>("BacktrackModuleLabel", "gaushitTruthMatch");
    _FlashMatchModuleLabel = pset.get<std::string>("FlashMatchModuleLabel", "flashmatch");
}

void SliceAnalysis::configure(fhicl::ParameterSet const &pset)
{
}

void SliceAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
    if (!e.getByLabel(_MCParticleModuleLabel, _mc_particle_handle))
        throw cet::exception("SliceAnalysis") << "failed to find any mc particles in event" << std::endl;
    art::fill_ptr_vector(_mc_particle_vector, _mc_particle_handle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(_mc_particle_vector, _mc_particle_map);

    if (!e.getByLabel(_PandoraModuleLabel, _pf_particle_handle))
        throw cet::exception("SliceAnalysis") << "failed to find any pandora-slice pf particles in event" << std::endl;
    art::fill_ptr_vector(_pf_particle_vector, _pf_particle_handle);
    lar_pandora::LArPandoraHelper::BuildPFParticleMap(_pf_particle_vector, _pf_particle_map);

    if (!e.getByLabel(_HitModuleLabel, _hit_handle))
        throw cet::exception("SliceAnalysis") << "failed ot find any hits in event" << std::endl;
    art::fill_ptr_vector(_hit_vector, _hit_handle);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(_hit_handle, e, _BacktrackModuleLabel);

    if (!e.getByLabel(_PandoraModuleLabel, _slice_handle))
        throw cet::exception("SliceAnalysis") << "failed to find any pandora slices in event" << std::endl;

    art::FindManyP<recob::Hit> hit_assoc = art::FindManyP<recob::Hit>(_slice_handle, e, _PandoraModuleLabel);
    art::fill_ptr_vector(_slice_vector, _slice_handle);

    art::FindManyP<recob::PFParticle> pf_part_slice_assoc = art::FindManyP<recob::PFParticle>(_slice_handle, e, _PandoraModuleLabel);
    art::FindManyP<larpandoraobj::PFParticleMetadata> pf_part_metadata_assoc = art::FindManyP<larpandoraobj::PFParticleMetadata>(_pf_particle_handle, e, _PandoraModuleLabel);

    if (!e.getByLabel(_FlashMatchModuleLabel, _flash_match_pfp_handle))
        throw cet::exception("SliceAnalysis") << "failed to find any falsh-matched pf particles" << std::endl;

    art::fill_ptr_vector(_flash_match_pfp_vector, _flash_match_pfp_handle);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(_flash_match_pfp_vector, _flash_nu_pfp_vector);

    // find mc particle information
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

            const int track_id = isParticleElectromagnetic(matched_mc_part) ? getLeadElectromagneticTrack(matched_mc_part) : matched_mc_part->TrackId();

            _hit_to_track_id[hit.key()] = track_id;
            _track_id_to_hits[track_id].push_back(hit);
        }
    }

    // find slice overall information and true neutrino slice
    int highest_hit_number(-1);
    std::map<int, int> slice_signal_hit_map;
    int total_true_hits(0);

    for (art::Ptr<recob::Slice> &slice : _slice_vector)
    {
        slice_signal_hit_map[slice->ID()] = 0; 

        const std::vector<art::Ptr<recob::Hit>> &slice_hit_vector(hit_assoc.at(slice.key()));
        for (const art::Ptr<recob::Hit> &slice_hit: slice_hit_vector)
        {
            if (_hit_to_track_id.find(slice_hit.key()) == _hit_to_track_id.end())
                continue;

            ++slice_signal_hit_map[slice->ID()];
            ++total_true_hits;
        }

        if ((slice_signal_hit_map[slice->ID()] > highest_hit_number) && (slice_signal_hit_map[slice->ID()] > 0))
        {
            highest_hit_number = slice_signal_hit_map[slice->ID()];
            _true_nu_slice_id = slice->ID();
        }
    }

    // find pandora slice 
    double best_topological_score(-std::numeric_limits<double>::max());

    for (const art::Ptr<recob::Slice> &slice: _slice_vector)
    {
        const std::vector<art::Ptr<recob::PFParticle>> pfp_slice_vector = pf_part_slice_assoc.at(slice.key());

        for (const art::Ptr<recob::PFParticle> &pfp : pfp_slice_vector)
        {
            if (!pfp->IsPrimary())
                continue;

            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfp_meta = pf_part_metadata_assoc.at(pfp.key());

            if (pfp_meta.empty())
                continue;

            const larpandoraobj::PFParticleMetadata::PropertiesMap &pf_particle_properties_map = pfp_meta.at(0)->GetPropertiesMap();
            if (!pf_particle_properties_map.empty() && (pf_particle_properties_map.find("NuScore") != pf_particle_properties_map.end()))
            {
                const double topological_score = pf_particle_properties_map.at("NuScore");
                
                if (topological_score > best_topological_score)
                {
                    best_topological_score = topological_score;
                    _pandora_nu_slice_id = slice->ID();
                }
            }
        }
    }

    // find flash match slice
    if (_flash_nu_pfp_vector.size() > 1) 
    {
        throw cet::exception("SliceAnalysis") << "Too many neutrinos found" << std::endl;
    }
    else if (_flash_nu_pfp_vector.size() == 1)
    {
        bool found = false;
    
        art::FindManyP<recob::Slice> flash_match_slice_assoc = art::FindManyP<recob::Slice>(_flash_match_pfp_handle, e, _FlashMatchModuleLabel);
        const std::vector<art::Ptr<recob::Slice>> &flash_match_slice_vector = flash_match_slice_assoc.at(_flash_nu_pfp_vector[0].key());

        if (!flash_match_slice_vector.empty())
        {
            const art::Ptr<recob::Slice> &flash_match_slice(flash_match_slice_vector.at(0));
            art::Handle<std::vector<recob::Slice>> flash_match_slice_handle;
            if (!e.getByLabel(_FlashMatchModuleLabel, flash_match_slice_handle))
                throw cet::exception("SliceAnalysis") << "failed to find flash match slice data" << std::endl;

            art::FindManyP<recob::Hit> flash_match_hit_assoc = art::FindManyP<recob::Hit>(_flash_match_pfp_handle, e, _FlashMatchModuleLabel);
            const std::vector<art::Ptr<recob::Hit>> &flash_match_slice_hits(flash_match_hit_assoc.at(flash_match_slice.key()));

            for (art::Ptr<recob::Slice> &slice : _slice_vector)
            {
                const std::vector<art::Ptr<recob::Hit>> &slice_hit_vector(hit_assoc.at(slice.key()));

                for (const art::Ptr<recob::Hit> &slice_hit : slice_hit_vector)
                {
                    for (const art::Ptr<recob::Hit> &flash_match_slice_hit : flash_match_slice_hits)
                    {
                        if (slice_hit.key() == flash_match_slice_hit.key())
                        {
                            found = true;
                            _flash_match_nu_slice_id = slice->ID();

                            break;
                        }
                    }

                    if (found)
                        break;
                }

                if (found)
                    break;
            }
        }
    }

    // for each slice type compute purity and completeness
    for (art::Ptr<recob::Slice> &slice : _slice_vector) 
    {
        const std::vector<art::Ptr<recob::Hit>> &slice_hits(hit_assoc.at(slice.key()));
        const int n_slice_hits = slice_hits.size();

        auto it = slice_signal_hit_map.find(slice->ID());
        const int n_slice_true_hits = (it == slice_signal_hit_map.end()) ? 0 : it->second;  

        float slice_completeness = total_true_hits == 0 ? 0.0 : static_cast<float>(n_slice_true_hits) / static_cast<float>(total_true_hits);
        float slice_purity = n_slice_hits == 0 ? 0.0 : static_cast<float>(n_slice_true_hits) / static_cast<float>(n_slice_hits);

        if (slice->ID() == _true_nu_slice_id && !fData)
        {
            _true_slice_completeness = slice_completeness;
            _true_slice_purity = slice_purity;
        }
        else if (slice->ID() == _pandora_nu_slice_id)
        {
            _pandora_slice_completeness = slice_completeness;
            _pandora_slice_purity = slice_purity;
        }
        else if (slice->ID() == _flash_match_nu_slice_id)
        {
            _flash_slice_completeness = slice_completeness;
            _flash_slice_purity = slice_purity;
        }
    }

    // find neutrino vertex information for each slice
    for (art::Ptr<recob::Slice> &slice : _slice_vector) 
    {
        const std::vector<art::Ptr<recob::PFParticle>> &slice_pf_part_vector(pf_part_slice_assoc.at(slice.key()));
        std::vector<art::Ptr<recob::PFParticle>> nu_pfp_vector;
        lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(slice_pf_part_vector, nu_pfp_vector);

        if (nu_pfp_vector.size() != 1)
            return;

        art::FindManyP<recob::Vertex> vertex_assoc = art::FindManyP<recob::Vertex>(_pf_particle_handle, e, _PandoraModuleLabel);
        const std::vector<art::Ptr<recob::Vertex>> &nu_vertex(vertex_assoc.at(nu_pfp_vector.at(0).key()));

        if (nu_vertex.empty())
            return;

        if (slice->ID() == _true_nu_slice_id && !fData)
        {
            _true_reco_nu_vtx_x = nu_vertex.at(0)->position().X();
            _true_reco_nu_vtx_y = nu_vertex.at(0)->position().Y();
            _true_reco_nu_vtx_z = nu_vertex.at(0)->position().Z();
        }
        else if (slice->ID() == _pandora_nu_slice_id)
        {
            _pandora_reco_nu_vtx_x = nu_vertex.at(0)->position().X();
            _pandora_reco_nu_vtx_y = nu_vertex.at(0)->position().Y();
            _pandora_reco_nu_vtx_z = nu_vertex.at(0)->position().Z();
        }
        else if (slice->ID() == _flash_match_nu_slice_id)
        {
            _flash_reco_nu_vtx_x = nu_vertex.at(0)->position().X();
            _flash_reco_nu_vtx_y = nu_vertex.at(0)->position().Y();
            _flash_reco_nu_vtx_z = nu_vertex.at(0)->position().Z();
        }
    }
}

void SliceAnalysis::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
}

void SliceAnalysis::setBranches(TTree *_tree)
{
    _tree->Branch("true_nu_slice_id", &_true_nu_slice_id, "true_nu_slice_id/I");
    _tree->Branch("pandora_nu_slice_id", &_pandora_nu_slice_id, "pandora_nu_slice_id/I");
    _tree->Branch("flash_match_nu_slice_id", &_flash_match_nu_slice_id, "flash_match_nu_slice_id/I");

    _tree->Branch("true_slice_completeness", &_true_slice_completeness, "true_slice_completeness/F");
    _tree->Branch("true_slice_purity", &_true_slice_purity, "true_slice_purity/F");
    _tree->Branch("pandora_slice_completeness", &_pandora_slice_completeness, "pandora_slice_completeness/F");
    _tree->Branch("pandora_slice_purity", &_pandora_slice_purity, "pandora_slice_purity/F");
    _tree->Branch("flash_slice_completeness", &_flash_slice_completeness, "flash_slice_completeness/F");
    _tree->Branch("flash_slice_purity", &_flash_slice_purity, "flash_slice_purity/F");

    _tree->Branch("true_reco_nu_vtx_x", &_true_reco_nu_vtx_x, "true_reco_nu_vtx_x/F");
    _tree->Branch("true_reco_nu_vtx_y", &_true_reco_nu_vtx_y, "true_reco_nu_vtx_y/F");
    _tree->Branch("true_reco_nu_vtx_z", &_true_reco_nu_vtx_z, "true_reco_nu_vtx_z/F");
    _tree->Branch("pandora_reco_nu_vtx_x", &_pandora_reco_nu_vtx_x, "pandora_reco_nu_vtx_x/F");
    _tree->Branch("pandora_reco_nu_vtx_y", &_pandora_reco_nu_vtx_y, "pandora_reco_nu_vtx_y/F");
    _tree->Branch("pandora_reco_nu_vtx_z", &_pandora_reco_nu_vtx_z, "pandora_reco_nu_vtx_z/F");
    _tree->Branch("flash_reco_nu_vtx_x", &_flash_reco_nu_vtx_x, "flash_reco_nu_vtx_x/F");
    _tree->Branch("flash_reco_nu_vtx_y", &_flash_reco_nu_vtx_y, "flash_reco_nu_vtx_y/F");
    _tree->Branch("flash_reco_nu_vtx_z", &_flash_reco_nu_vtx_z, "flash_reco_nu_vtx_z/F");
}

void SliceAnalysis::resetTTree(TTree *_tree)
{
    _true_nu_slice_id = -1;
    _pandora_nu_slice_id = -1;
    _flash_match_nu_slice_id = -1;

    _true_slice_completeness = std::numeric_limits<float>::lowest();
    _true_slice_purity = std::numeric_limits<float>::lowest();
    _pandora_slice_completeness = std::numeric_limits<float>::lowest();
    _pandora_slice_purity = std::numeric_limits<float>::lowest();
    _flash_slice_completeness = std::numeric_limits<float>::lowest();
    _flash_slice_purity = std::numeric_limits<float>::lowest();

    _true_reco_nu_vtx_x = std::numeric_limits<float>::lowest();
    _true_reco_nu_vtx_y = std::numeric_limits<float>::lowest();
    _true_reco_nu_vtx_z = std::numeric_limits<float>::lowest();
    _pandora_reco_nu_vtx_x = std::numeric_limits<float>::lowest();
    _pandora_reco_nu_vtx_y = std::numeric_limits<float>::lowest();
    _pandora_reco_nu_vtx_z = std::numeric_limits<float>::lowest();
    _flash_reco_nu_vtx_x = std::numeric_limits<float>::lowest();
    _flash_reco_nu_vtx_y = std::numeric_limits<float>::lowest();
    _flash_reco_nu_vtx_z = std::numeric_limits<float>::lowest();
}

bool SliceAnalysis::isParticleElectromagnetic(const art::Ptr<simb::MCParticle> &mc_part)
{
    return ((std::abs(mc_part->PdgCode() == 11) || (mc_part->PdgCode() == 22)));
}

int SliceAnalysis::getLeadElectromagneticTrack(const art::Ptr<simb::MCParticle> &mc_part)
{
    int track_id = mc_part->TrackId();
    art::Ptr<simb::MCParticle> mother_mc_part = mc_part;

    do 
    {
        track_id = mother_mc_part->TrackId();
        const int mother_id = mother_mc_part->Mother();

        if (_mc_particle_map.find(mother_id) == _mc_particle_map.end())
            break;

        mother_mc_part = _mc_particle_map.at(mother_id);
    } while (isParticleElectromagnetic(mother_mc_part));

    return track_id;
}

DEFINE_ART_CLASS_TOOL(SliceAnalysis)
} 

#endif