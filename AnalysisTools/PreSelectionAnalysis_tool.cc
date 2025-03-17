#ifndef ANALYSIS_PRESELECTION_CXX
#define ANALYSIS_PRESELECTION_CXX

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
#include "lardataobj/RecoBase/OpFlash.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

namespace analysis
{

class PreSelectionAnalysis : public AnalysisToolBase
{

public:
    PreSelectionAnalysis(const fhicl::ParameterSet &pset);
    ~PreSelectionAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);
    void analyseEvent(art::Event const &e, bool fData) override;
    void analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;
    void SaveTruth(art::Event const &e);
    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

private:
    art::InputTag _SimulationModuleLabel;
    art::InputTag _PandoraModuleLabel;
    art::InputTag _HitModuleLabel;
    art::InputTag _BacktrackModuleLabel;
    art::InputTag _FlashMatchModuleLabel;
    art::InputTag _SpacePointModuleLabel;
    art::InputTag _PFParticleModuleLabel;
    art::InputTag _FlashLabel;

    TTree* _slice_tree;

    std::vector<float> _flash_time_v;
    std::vector<float> _flash_total_pe_v;
    std::vector<float> _flash_center_y_v;
    std::vector<float> _flash_center_z_v;
    std::vector<float> _flash_width_y_v;
    std::vector<float> _flash_width_z_v;

    std::vector<float> _slice_completeness_v;
    std::vector<float> _slice_purity_v;
    std::vector<float> _slice_topological_score_v;
    std::vector<float> _slice_pandora_score_v;
    std::vector<int> _slice_ids_v;
    std::vector<int> _slice_n_hits_v;

    std::vector<float> _slice_center_x_v;
    std::vector<float> _slice_center_y_v;
    std::vector<float> _slice_center_z_v;
    std::vector<float> _slice_charge_v;

    int _true_nu_slice_id;
    int _pandora_nu_slice_id;
    int _flash_match_nu_slice_id;

    bool _flash_slice_found;
    bool _pandora_slice_found;

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

    bool _debug;
};

PreSelectionAnalysis::PreSelectionAnalysis(const fhicl::ParameterSet &pset)
    : _debug(pset.get<bool>("DebugMode", false))
{
    _SimulationModuleLabel = pset.get<std::string>("MCParticleModuleLabel", "largeant");
    _PandoraModuleLabel = pset.get<std::string>("PandoraModuleLabel", "pandora");
    _HitModuleLabel = pset.get<std::string>("HitModuleLabel", "gaushit");
    _BacktrackModuleLabel = pset.get<std::string>("BacktrackModuleLabel", "gaushitTruthMatch");
    _FlashMatchModuleLabel = pset.get<std::string>("FlashMatchModuleLabel", "pandora");
    _PFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel", "pandora");
    _SpacePointModuleLabel = pset.get<std::string>("SpacePointModuleLabel", "pandora"); 
    _FlashLabel = pset.get<std::string>("FlashLabel", "simpleFlashBeam");

    art::ServiceHandle<art::TFileService> tfs;
    _slice_tree = tfs->make<TTree>("PreSelectionAnalysis", "Slice Tree");
}

void PreSelectionAnalysis::configure(fhicl::ParameterSet const &pset)
{
}

void PreSelectionAnalysis::analyseEvent(art::Event const &e, bool fData)
{
    std::cout << "Analysing event in PreSelection..." << std::endl;
    art::Handle<std::vector<simb::MCParticle>> mc_particle_handle;
    std::vector<art::Ptr<simb::MCParticle>> mc_particle_vector;

    lar_pandora::PFParticleMap pf_particle_map;
    lar_pandora::MCParticleMap mc_particle_map; 

    if (!e.getByLabel(_SimulationModuleLabel, mc_particle_handle))
        throw cet::exception("PreSelectionAnalysis") << "Failed to find any MC particles in event" << std::endl;
    art::fill_ptr_vector(mc_particle_vector, mc_particle_handle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mc_particle_vector, mc_particle_map);

    art::Handle<std::vector<recob::PFParticle>> pf_particle_handle;
    std::vector<art::Ptr<recob::PFParticle>> pf_particle_vector;

    if (!e.getByLabel(_PandoraModuleLabel, pf_particle_handle))
        throw cet::exception("PreSelectionAnalysis") << "Failed to find any Pandora-slice PFParticles in event" << std::endl;
    art::fill_ptr_vector(pf_particle_vector, pf_particle_handle);
    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pf_particle_vector, pf_particle_map);

    art::Handle<std::vector<recob::Hit>> hit_handle;
    std::vector<art::Ptr<recob::Hit>> hit_vector;

    if (!e.getByLabel(_HitModuleLabel, hit_handle))
        throw cet::exception("PreSelectionAnalysis") << "Failed to find any hits in event" << std::endl;
    art::fill_ptr_vector(hit_vector, hit_handle);

    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hit_handle, e, _BacktrackModuleLabel);

    art::Handle<std::vector<recob::Slice>> slice_handle; 
    std::vector<art::Ptr<recob::Slice>> slice_vector;
    if (!e.getByLabel(_PandoraModuleLabel, slice_handle))
        throw cet::exception("PreSelectionAnalysis") << "Failed to find any Pandora slices in event" << std::endl;

    art::FindManyP<recob::Hit> hit_slice_assoc = art::FindManyP<recob::Hit>(slice_handle, e, _PandoraModuleLabel);
    art::FindManyP<recob::SpacePoint> sp_slice_assoc = art::FindManyP<recob::SpacePoint>(slice_handle, e, _PandoraModuleLabel);
    art::fill_ptr_vector(slice_vector, slice_handle);

    art::FindManyP<recob::PFParticle> pf_part_slice_assoc = art::FindManyP<recob::PFParticle>(slice_handle, e, _PandoraModuleLabel);
    art::FindManyP<larpandoraobj::PFParticleMetadata> pf_part_metadata_assoc = art::FindManyP<larpandoraobj::PFParticleMetadata>(pf_particle_handle, e, _PandoraModuleLabel);

    art::Handle<std::vector<recob::PFParticle>> flash_match_pfp_handle;
    std::vector<art::Ptr<recob::PFParticle>> flash_match_pfp_vector;
    std::vector<art::Ptr<recob::PFParticle>> flash_nu_pfp_vector;

    if (!e.getByLabel(_FlashMatchModuleLabel, flash_match_pfp_handle))
        throw cet::exception("PreSelectionAnalysis") << "Failed to find any flash-matched PFParticles" << std::endl;

    art::fill_ptr_vector(flash_match_pfp_vector, flash_match_pfp_handle);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(flash_match_pfp_vector, flash_nu_pfp_vector);

    art::InputTag flash_tag(_FlashLabel);
    const auto flashes(*e.getValidHandle<std::vector<recob::OpFlash>>(flash_tag));

    auto const& pfp_handle = e.getValidHandle<std::vector<recob::PFParticle>>(_PFParticleModuleLabel);
    auto const& sp_handle = e.getValidHandle<std::vector<recob::SpacePoint>>(_SpacePointModuleLabel);
    art::FindManyP<recob::SpacePoint> sp_pfp_assoc(pfp_handle, e, _PFParticleModuleLabel);
    art::FindManyP<recob::Hit> hit_sp_assoc(sp_handle, e, _SpacePointModuleLabel);

    for (const auto& op_flash : flashes)
    {
        _flash_time_v.push_back(op_flash.Time());
        _flash_total_pe_v.push_back(op_flash.TotalPE());
        _flash_center_y_v.push_back(op_flash.YCenter());
        _flash_center_z_v.push_back(op_flash.ZCenter());
        _flash_width_y_v.push_back(op_flash.YWidth());
        _flash_width_z_v.push_back(op_flash.ZWidth());
    }

    std::map<int, int> hit_to_track_id;
    std::map<int, std::vector<art::Ptr<recob::Hit>>> track_id_to_hits;

    // find simulation particle information
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

            hit_to_track_id[hit.key()] = track_idx;
            track_id_to_hits[track_idx].push_back(hit);
        }
    }

    // find slice overall information and true neutrino slice
    int highest_hit_number(-1);
    std::map<int, int> slice_signal_hit_map;
    int total_true_hits(0);

    for (art::Ptr<recob::Slice> &slice : slice_vector)
    {
        slice_signal_hit_map[slice->ID()] = 0; 

        const std::vector<art::Ptr<recob::Hit>> &slice_hit_vector(hit_slice_assoc.at(slice.key()));

        for (const art::Ptr<recob::Hit> &slice_hit : slice_hit_vector)
        {
            if (hit_to_track_id.find(slice_hit.key()) == hit_to_track_id.end()) 
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

    // find Pandora slice 
    double best_topological_score(-std::numeric_limits<double>::max());

    for (const art::Ptr<recob::Slice> &slice: slice_vector)
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
                    _pandora_slice_found = true;
                }
            }
        }
    }

    // find flash match slice
    if (flash_nu_pfp_vector.size() > 1) 
    {
        throw cet::exception("PreSelectionAnalysis") << "Too many neutrinos found" << std::endl;
    }
    else if (flash_nu_pfp_vector.size() == 1)
    {
        //bool found = false;

        art::FindManyP<recob::Slice> flash_match_slice_assoc = art::FindManyP<recob::Slice>(flash_match_pfp_handle, e, _FlashMatchModuleLabel);
        const std::vector<art::Ptr<recob::Slice>> &flash_match_slice_vector = flash_match_slice_assoc.at(flash_nu_pfp_vector[0].key());

        if (!flash_match_slice_vector.empty())
        {
            const art::Ptr<recob::Slice> &flash_match_slice = flash_match_slice_vector[0];
            _flash_match_nu_slice_id = flash_match_slice->ID();
            _flash_slice_found = true;
        } 
    }

    // For each slice type compute purity and completeness
    for (art::Ptr<recob::Slice> &slice : slice_vector) 
    {
        const std::vector<art::Ptr<recob::Hit>> &slice_hits(hit_slice_assoc.at(slice.key()));
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
    for (art::Ptr<recob::Slice> &slice : slice_vector) 
    {
        const std::vector<art::Ptr<recob::PFParticle>> &slice_pf_part_vector(pf_part_slice_assoc.at(slice.key()));
        std::vector<art::Ptr<recob::PFParticle>> nu_pfp_vector;
        lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(slice_pf_part_vector, nu_pfp_vector);

        if (nu_pfp_vector.size() != 1)
            continue;

        art::FindManyP<recob::Vertex> vertex_assoc = art::FindManyP<recob::Vertex>(pf_particle_handle, e, _PandoraModuleLabel);
        const std::vector<art::Ptr<recob::Vertex>> &nu_vertex(vertex_assoc.at(nu_pfp_vector.at(0).key()));

        if (nu_vertex.empty())
            continue;

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

    // for each slice, find each information
    for (art::Ptr<recob::Slice> &slice : slice_vector)
    {
        const std::vector<art::Ptr<recob::Hit>> &slice_hits(hit_slice_assoc.at(slice.key()));
        const int n_slice_hits = slice_hits.size();

        auto it = slice_signal_hit_map.find(slice->ID());
        const int n_slice_true_hits = (it == slice_signal_hit_map.end()) ? 0 : it->second;

        float slice_completeness = total_true_hits == 0 ? 0.0 : static_cast<float>(n_slice_true_hits) / static_cast<float>(total_true_hits);
        float slice_purity = n_slice_hits == 0 ? 0.0 : static_cast<float>(n_slice_true_hits) / static_cast<float>(n_slice_hits);

        _slice_completeness_v.push_back(slice_completeness);
        _slice_purity_v.push_back(slice_purity);
        _slice_ids_v.push_back(slice->ID());
        _slice_n_hits_v.push_back(n_slice_hits);

        double topological_score = -std::numeric_limits<double>::max();
        double pandora_score = -std::numeric_limits<double>::max();

        const std::vector<art::Ptr<recob::PFParticle>> pfp_slice_vector = pf_part_slice_assoc.at(slice.key());

        for (const art::Ptr<recob::PFParticle> &pfp : pfp_slice_vector)
        {
            if (!pfp->IsPrimary())
                continue;

            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfp_meta = pf_part_metadata_assoc.at(pfp.key());

            if (!pfp_meta.empty())
            {
                const larpandoraobj::PFParticleMetadata::PropertiesMap &pf_particle_properties_map = pfp_meta.at(0)->GetPropertiesMap();
                
                if (pf_particle_properties_map.find("NuScore") != pf_particle_properties_map.end())
                {
                    topological_score = pf_particle_properties_map.at("NuScore");
                }

                if (pf_particle_properties_map.find("PandoraScore") != pf_particle_properties_map.end())
                {
                    pandora_score = pf_particle_properties_map.at("PandoraScore");
                }
            }
        }

        _slice_topological_score_v.push_back(topological_score);
        _slice_pandora_score_v.push_back(pandora_score);

        std::vector<art::Ptr<recob::SpacePoint>> slice_spacepoints_v;
        std::vector<std::vector<art::Ptr<recob::Hit>>> slice_hits_v;

        for (const art::Ptr<recob::PFParticle> &pfp : pfp_slice_vector)
        {
            if (pfp->IsPrimary())
                continue;

            try
            {
                const std::vector<art::Ptr<recob::SpacePoint>>& pfp_spacepoints_v = sp_pfp_assoc.at(pfp.key());

                std::vector<art::Ptr<recob::Hit>> pfp_hits_v;
                for (const art::Ptr<recob::SpacePoint> sp : pfp_spacepoints_v)
                {
                    try
                    {
                        const std::vector<art::Ptr<recob::Hit>> sp_hits_v = hit_sp_assoc.at(sp.key());

                        slice_spacepoints_v.push_back(sp);
                        slice_hits_v.push_back(sp_hits_v); 
                    }
                    catch (const std::out_of_range& e)
                    {
                        continue; 
                    }
                }
            }
            catch (const std::out_of_range& e)
            {
                continue; 
            }
        }

        float total_charge = 0.0;
        float weighted_x = 0.0;
        float weighted_y = 0.0;
        float weighted_z = 0.0;

        for (size_t i = 0; i < slice_spacepoints_v.size(); ++i)
        {
            const art::Ptr<recob::SpacePoint> &spacepoint = slice_spacepoints_v[i];
            const std::vector<art::Ptr<recob::Hit>> &hit_v = slice_hits_v[i];

            float charge = 0;

            for (const art::Ptr<recob::Hit> hit : hit_v)
            {
                if (hit->View() != geo::kZ)
	                continue;

                charge += hit->Integral();
            }

            const double *xyz = spacepoint->XYZ();
            total_charge += charge;

            weighted_x += charge * xyz[0];
            weighted_y += charge * xyz[1];
            weighted_z += charge * xyz[2];
        }

        if (total_charge > 0.0)
        {
            float center_x = weighted_x / total_charge;
            float center_y = weighted_y / total_charge;
            float center_z = weighted_z / total_charge;

            _slice_center_x_v.push_back(center_x);
            _slice_center_y_v.push_back(center_y);
            _slice_center_z_v.push_back(center_z);
        }
        else
        {
            _slice_center_x_v.push_back(0.0);
            _slice_center_y_v.push_back(0.0);
            _slice_center_z_v.push_back(0.0);
        }

        _slice_charge_v.push_back(total_charge);        
    }

    _slice_tree->Fill();

    std::cout << "Finished analysing event in PreSelection!" << std::endl;
}

void PreSelectionAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
}

void PreSelectionAnalysis::setBranches(TTree *_tree)
{
    _slice_tree->Branch("flash_time", &_flash_time_v);
    _slice_tree->Branch("flash_total_pe", &_flash_total_pe_v);
    _slice_tree->Branch("flash_center_y", &_flash_center_y_v);
    _slice_tree->Branch("flash_center_z", &_flash_center_z_v);
    _slice_tree->Branch("flash_width_y", &_flash_width_y_v);
    _slice_tree->Branch("flash_width_z", &_flash_width_z_v);

    _slice_tree->Branch("slice_completeness", &_slice_completeness_v);
    _slice_tree->Branch("slice_purity", &_slice_purity_v);
    _slice_tree->Branch("slice_topological_score", &_slice_topological_score_v);
    _slice_tree->Branch("slice_pandora_score", &_slice_pandora_score_v);
    _slice_tree->Branch("slice_ids", &_slice_ids_v);
    _slice_tree->Branch("slice_n_hits", &_slice_n_hits_v);

    _slice_tree->Branch("slice_center_x", &_slice_center_x_v);
    _slice_tree->Branch("slice_center_y", &_slice_center_y_v);
    _slice_tree->Branch("slice_center_z", &_slice_center_z_v);
    _slice_tree->Branch("slice_charge", &_slice_charge_v);

    _slice_tree->Branch("true_nu_slice_id", &_true_nu_slice_id, "true_nu_slice_id/I");
    _slice_tree->Branch("pandora_nu_slice_id", &_pandora_nu_slice_id, "pandora_nu_slice_id/I");
    _slice_tree->Branch("flash_match_nu_slice_id", &_flash_match_nu_slice_id, "flash_match_nu_slice_id/I");

    _slice_tree->Branch("true_slice_completeness", &_true_slice_completeness, "true_slice_completeness/F");
    _slice_tree->Branch("true_slice_purity", &_true_slice_purity, "true_slice_purity/F");
    _slice_tree->Branch("pandora_slice_completeness", &_pandora_slice_completeness, "pandora_slice_completeness/F");
    _slice_tree->Branch("pandora_slice_purity", &_pandora_slice_purity, "pandora_slice_purity/F");
    _slice_tree->Branch("flash_slice_completeness", &_flash_slice_completeness, "flash_slice_completeness/F");
    _slice_tree->Branch("flash_slice_purity", &_flash_slice_purity, "flash_slice_purity/F");

    _slice_tree->Branch("true_reco_nu_vtx_x", &_true_reco_nu_vtx_x, "true_reco_nu_vtx_x/F");
    _slice_tree->Branch("true_reco_nu_vtx_y", &_true_reco_nu_vtx_y, "true_reco_nu_vtx_y/F");
    _slice_tree->Branch("true_reco_nu_vtx_z", &_true_reco_nu_vtx_z, "true_reco_nu_vtx_z/F");
    _slice_tree->Branch("pandora_reco_nu_vtx_x", &_pandora_reco_nu_vtx_x, "pandora_reco_nu_vtx_x/F");
    _slice_tree->Branch("pandora_reco_nu_vtx_y", &_pandora_reco_nu_vtx_y, "pandora_reco_nu_vtx_y/F");
    _slice_tree->Branch("pandora_reco_nu_vtx_z", &_pandora_reco_nu_vtx_z, "pandora_reco_nu_vtx_z/F");
    _slice_tree->Branch("flash_reco_nu_vtx_x", &_flash_reco_nu_vtx_x, "flash_reco_nu_vtx_x/F");
    _slice_tree->Branch("flash_reco_nu_vtx_y", &_flash_reco_nu_vtx_y, "flash_reco_nu_vtx_y/F");
    _slice_tree->Branch("flash_reco_nu_vtx_z", &_flash_reco_nu_vtx_z, "flash_reco_nu_vtx_z/F");

    _slice_tree->Branch("flash_slice_found", &_flash_slice_found, "flash_slice_found/O"); 
    _slice_tree->Branch("pandora_slice_found", &_pandora_slice_found, "pandora_slice_found/O");
}

void PreSelectionAnalysis::resetTTree(TTree *_tree)
{
    _flash_time_v.clear();
    _flash_total_pe_v.clear();
    _flash_center_y_v.clear();
    _flash_center_z_v.clear();
    _flash_width_y_v.clear();
    _flash_width_z_v.clear();

    _slice_completeness_v.clear();
    _slice_purity_v.clear();
    _slice_topological_score_v.clear();
    _slice_pandora_score_v.clear();
    _slice_ids_v.clear();
    _slice_n_hits_v.clear();

    _slice_center_x_v.clear();
    _slice_center_y_v.clear();
    _slice_center_z_v.clear();
    _slice_charge_v.clear();

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

    _flash_slice_found = false;
    _pandora_slice_found = false;
}

DEFINE_ART_CLASS_TOOL(PreSelectionAnalysis)
} 

#endif