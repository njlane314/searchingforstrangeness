#ifndef DEFAULT_ANALYSIS_CXX
#define DEFAULT_ANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "Common/BacktrackingFuncs.h"
#include "Common/Geometry.h"
#include "Common/SCECorrections.h"
#include "Common/Containment.h"
#include "Common/TrackShowerScoreFuncs.h"
#include "Common/ProximityClustering.h"
#include "Common/Descendents.h"
#include "Common/Scatters.h"
#include "Common/Pandora.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RawData/TriggerData.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"

namespace analysis {

class DefaultAnalysis : public AnalysisToolBase {
public:
    DefaultAnalysis(const fhicl::ParameterSet &pset);
    ~DefaultAnalysis(){};
    void configure(fhicl::ParameterSet const &pset);
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool is_data, bool is_selected) override;
    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

private:
    art::InputTag fPFPproducer; 
    art::InputTag fCRTproducer;
    art::InputTag fCLSproducer;
    art::InputTag fMCTproducer;
    art::InputTag fMCPproducer;
    art::InputTag fMCFproducer;
    art::InputTag fBKTproducer;
    art::InputTag fHITproducer;
    art::InputTag fMCRproducer;
    art::InputTag fSLCproducer;

    float fTrkShrScore;
    std::string NuMIOpFilterProd;
    std::string NuMISWTrigProd;

    float fFidvolXstart;
    float fFidvolXend;
    float fFidvolYstart;
    float fFidvolYend;
    float fFidvolZstart;
    float fFidvolZend;

    bool fMakeNuMINtuple;
    bool fIgnoreMCFlux;

    std::vector<float> _slice_topological_scores;
    float _reco_neutrino_vertex_x;
    float _reco_neutrino_vertex_y;
    float _reco_neutrino_vertex_z;
    float _reco_neutrino_vertex_sce_x;
    float _reco_neutrino_vertex_sce_y;
    float _reco_neutrino_vertex_sce_z;

    int _software_trigger;
    int _software_trigger_pre;
    int _software_trigger_post;
    int _software_trigger_pre_ext;
    int _software_trigger_post_ext;

    float _optical_filter_pe_beam;
    float _optical_filter_pe_veto;

    int _num_slices;

    int _crt_veto;
    float _crt_hit_pe;

    std::vector<int> _pfp_slice_indices;
    std::vector<int> _backtracked_pdg_codes;
    std::vector<float> _backtracked_energies;
    std::vector<int> _backtracked_track_ids;
    std::vector<float> _backtracked_purities;
    std::vector<float> _backtracked_completenesses;
    std::vector<float> _backtracked_overlay_purities;
    std::vector<float> _backtracked_momentum_x;
    std::vector<float> _backtracked_momentum_y;
    std::vector<float> _backtracked_momentum_z;
    std::vector<float> _backtracked_start_x;
    std::vector<float> _backtracked_start_y;
    std::vector<float> _backtracked_start_z;
    std::vector<float> _backtracked_start_time;
    std::vector<float> _backtracked_start_wire_U;
    std::vector<float> _backtracked_start_wire_V;
    std::vector<float> _backtracked_start_wire_Y;
    std::vector<float> _backtracked_sce_start_x;
    std::vector<float> _backtracked_sce_start_y;
    std::vector<float> _backtracked_sce_start_z;
    std::vector<float> _backtracked_sce_start_wire_U;
    std::vector<float> _backtracked_sce_start_wire_V;
    std::vector<float> _backtracked_sce_start_wire_Y;

    int _selection_pass;
    int _event_total_hits;
    int _slice_pdg;
    int _slice_num_hits;
    int _slice_id;
    float _topological_score;

    std::vector<int> _pfp_pdg_codes;
    std::vector<int> _pfp_num_hits;
    std::vector<int> _pfp_num_plane_hits_U;
    std::vector<int> _pfp_num_plane_hits_V;
    std::vector<int> _pfp_num_plane_hits_Y;
    std::vector<int> _pfp_num_subclusters_U;
    std::vector<int> _pfp_num_subclusters_V;
    std::vector<int> _pfp_num_subclusters_Y;
    std::vector<float> _pfp_max_subhit_fraction_U;
    std::vector<float> _pfp_max_subhit_fraction_V;
    std::vector<float> _pfp_max_subhit_fraction_Y;

    float _slice_cluster_fraction;
    std::vector<uint> _pfp_generations;
    std::vector<uint> _pfp_shower_daughters;
    std::vector<uint> _pfp_track_daughters;
    std::vector<uint> _pfp_num_descendents;
    std::vector<float> _pfp_vertex_x;
    std::vector<float> _pfp_vertex_y;
    std::vector<float> _pfp_vertex_z;
    unsigned int _num_pfps;
    std::vector<float> _track_shower_scores;
    unsigned int _num_tracks;
    unsigned int _num_showers;
    unsigned int _total_hits_U;
    unsigned int _total_hits_V;
    unsigned int _total_hits_Y;
};

DefaultAnalysis::DefaultAnalysis(const fhicl::ParameterSet &p) {
    fPFPproducer = p.get<art::InputTag>("PFPproducer", "pandora");
    fCRTproducer = p.get<art::InputTag>("CRTproducer", "pandora");
    fCLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
    fMCTproducer = p.get<art::InputTag>("MCTproducer", "generator");
    fMCPproducer = p.get<art::InputTag>("MCPproducer", "largeant");
    fMCFproducer = p.get<art::InputTag>("MCFproducer", "generator");
    fBKTproducer = p.get<art::InputTag>("BKTproducer", "gaushitTruthMatch");
    fHITproducer = p.get<art::InputTag>("HITproducer", "gaushit");
    fMCRproducer = p.get<art::InputTag>("MCRproducer", "mcreco");
    fSLCproducer = p.get<art::InputTag>("SLCproducer", "pandora");

    fTrkShrScore = p.get<float>("TrkShrScore", 0.5);
    
    fFidvolXstart = p.get<double>("fidvolXstart", 10);
    fFidvolXend = p.get<double>("fidvolXend", 10);
    fFidvolYstart = p.get<double>("fidvolYstart", 15);
    fFidvolYend = p.get<double>("fidvolYend", 15);
    fFidvolZstart = p.get<double>("fidvolZstart", 10);
    fFidvolZend = p.get<double>("fidvolZend", 50);

    fMakeNuMINtuple = p.get<bool>("makeNuMINtuple", false);
    fIgnoreMCFlux = p.get<bool>("ignoreMCFlux", false);
    NuMIOpFilterProd = p.get<std::string>("NuMIOpFiltProcName", "DataStage1Optical");
    NuMISWTrigProd = p.get<std::string>("NuMISWTriggerProcName", "DataOverlayOpticalNuMI");
}

void DefaultAnalysis::configure(fhicl::ParameterSet const &p) {
}

void DefaultAnalysis::analyseEvent(const art::Event &event, bool is_data) {
    common::ProxySliceColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(event, fPFPproducer,
        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
        proxy::withAssociated<recob::Slice>(fSLCproducer));
    int pfp_slice_id;
    int temp_pfp_slice_id;
    int max_slice_id = 0;
    for (const common::ProxySliceElem_t &pfp : pfp_proxy) {
        auto temp_slice_pxy_v = pfp.get<recob::Slice>();
        if (temp_slice_pxy_v.size() != 0) {
            temp_pfp_slice_id = temp_slice_pxy_v.at(0)->ID();
            if (temp_pfp_slice_id > max_slice_id) {
                max_slice_id = temp_pfp_slice_id;
            }
        }
    }

    std::vector<float> temp_slice_topo_score_v(max_slice_id + 1);
    fill(temp_slice_topo_score_v.begin(), temp_slice_topo_score_v.end(), std::numeric_limits<float>::lowest());
    for (const common::ProxySliceElem_t &pfp : pfp_proxy) {
        auto metadata_pxy_v = pfp.get<larpandoraobj::PFParticleMetadata>();
        auto slice_pxy_v = pfp.get<recob::Slice>();
        if (slice_pxy_v.size() != 0) {
            pfp_slice_id = slice_pxy_v.at(0)->ID();
            if (metadata_pxy_v.size() != 0) {
                for (unsigned int j = 0; j < metadata_pxy_v.size(); ++j) {
                    const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(metadata_pxy_v.at(j));
                    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
                    if (!pfParticlePropertiesMap.empty() && temp_slice_topo_score_v.at(pfp_slice_id) == std::numeric_limits<float>::lowest()) {
                        auto it = pfParticlePropertiesMap.begin();
                        while (it != pfParticlePropertiesMap.end()) {
                            if (it->first == "NuScore") {
                                temp_slice_topo_score_v.at(pfp_slice_id) = pfParticlePropertiesMap.at(it->first);
                            }
                            it++;
                        }
                    }
                }
            }
        }
    }

    _slice_topological_scores = temp_slice_topo_score_v;

    if (!is_data && (!fMakeNuMINtuple)) {
        art::Handle<uboone::UbooneOpticalFilter> CommonOpticalFilter_h;
        art::InputTag fCommonOppFiltTag("opfiltercommon");
        event.getByLabel(fCommonOppFiltTag, CommonOpticalFilter_h);
        _optical_filter_pe_beam = CommonOpticalFilter_h->PE_Beam();
        _optical_filter_pe_veto = CommonOpticalFilter_h->PE_Veto();
    }

    if (fMakeNuMINtuple) {
        art::Handle<uboone::UbooneOpticalFilter> CommonOpticalFilter_h;
        std::cout << "Test1: " << NuMIOpFilterProd << std::endl;
        art::InputTag fCommonOppFiltTag("opfiltercommon", "", NuMIOpFilterProd);
        event.getByLabel(fCommonOppFiltTag, CommonOpticalFilter_h);
        if (!CommonOpticalFilter_h.isValid()) {
            std::cout << "Could not find data override for common op filter using default... (this is expected for overlay)" << std::endl;
            art::InputTag fCommonOppFiltTag("opfiltercommon");
            event.getByLabel(fCommonOppFiltTag, CommonOpticalFilter_h);
        }
        _optical_filter_pe_beam = CommonOpticalFilter_h->PE_Beam();
        _optical_filter_pe_veto = CommonOpticalFilter_h->PE_Veto();

        std::cout << "Test2: " << NuMISWTrigProd << std::endl;
        art::InputTag swtrig_tag("TriggerResults", "", NuMISWTrigProd);
        art::Handle<art::TriggerResults> swtrig_handle;
        event.getByLabel(swtrig_tag, swtrig_handle);
        if (swtrig_handle.isValid()) {
            if (swtrig_handle->accept() == true) {
                _software_trigger = 1;
            }
            else {
                _software_trigger = 0;
            }
        }

        if (!is_data) {
            art::InputTag triggerTag("swtrigger", "", NuMISWTrigProd);
            const auto& triggerHandle = event.getValidHandle<raw::ubdaqSoftwareTriggerData>(triggerTag);
            std::vector<std::string> triggerName = triggerHandle->getListOfAlgorithms();
            for (int j = 0; j != triggerHandle->getNumberOfAlgorithms(); j++) {
                if (triggerName[j] == "EXT_NUMIwin_FEMBeamTriggerAlgo") {
                    _software_trigger_pre_ext = triggerHandle->passedAlgo(triggerName[j]);
                }
                else if (triggerName[j] == "EXT_NUMIwin_2018May_FEMBeamTriggerAlgo") {
                    _software_trigger_post_ext = triggerHandle->passedAlgo(triggerName[j]);
                }
                else if (triggerName[j] == "NUMI_FEMBeamTriggerAlgo") {
                    _software_trigger_pre = triggerHandle->passedAlgo(triggerName[j]);
                }
                else if (triggerName[j] == "NUMI_2018May_FEMBeamTriggerAlgo") {
                    _software_trigger_post = triggerHandle->passedAlgo(triggerName[j]);
                }
                else {
                    continue;
                }
                std::cout << triggerName[j] << ": ";
                std::cout << triggerHandle->passedAlgo(triggerName[j]) << std::endl;
            }
        }
    }
    else {
        art::InputTag swtrig_tag("TriggerResults", "", "DataOverlayOptical");
        art::Handle<art::TriggerResults> swtrig_handle;
        event.getByLabel(swtrig_tag, swtrig_handle);
        if (swtrig_handle.isValid()) {
            if (swtrig_handle->accept() == true) {
                _software_trigger = 1;
            }
            else {
                _software_trigger = 0;
            }
        }
    }

    if (fCRTproducer != "") {
        art::Handle<art::Assns<crt::CRTHit, recob::OpFlash, void>> crtveto_h;
        event.getByLabel(fCRTproducer, crtveto_h);
        _crt_veto = crtveto_h->size();
        if (_crt_veto == 1) {
            _crt_hit_pe = crtveto_h->at(0).first->peshit;
        }
    }

    art::ValidHandle<std::vector<recob::Hit>> inputHits = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    _event_total_hits = inputHits->size();
}

void DefaultAnalysis::analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool is_data, bool is_selected) {
    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(event, fCLSproducer,
        proxy::withAssociated<recob::Hit>(fCLSproducer));
    art::ValidHandle<std::vector<recob::Slice>> inputSlice = event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, event, fSLCproducer));

    lar_pandora::LArPandoraHelper larpandora;
    lar_pandora::PFParticleVector pfparticles;
    lar_pandora::PFParticleMap particleMap;
    larpandora.CollectPFParticles(event, "pandora", pfparticles);
    larpandora.BuildPFParticleMap(pfparticles, particleMap);

    std::vector<common::BtPart> btparts_v;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
    if (!is_data) {
        const std::vector<sim::MCShower> &inputMCShower = *(event.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
        const std::vector<sim::MCTrack> &inputMCTrack = *(event.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
        art::ValidHandle<std::vector<recob::Hit>> inputHits = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
        assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, event, fBKTproducer));
        btparts_v = common::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);
    }

    size_t pfpidx = 0;
    _num_pfps = 0;
    for (auto pfp : slice_pfp_vec) {
        if (pfp->IsPrimary()) {
            _slice_pdg = pfp->PdgCode();
            auto slice_pxy_v = pfp.get<recob::Slice>();
            if (slice_pxy_v.size() != 1) 
                return;
            
            auto slicehits = assocSliceHit->at(slice_pxy_v[0].key());
            _slice_num_hits = slicehits.size();
            auto metadata_pxy_v = pfp.get<larpandoraobj::PFParticleMetadata>();
            _slice_id = slice_pxy_v.at(0)->ID();
            if (metadata_pxy_v.size() != 0) {
                for (unsigned int j = 0; j < metadata_pxy_v.size(); ++j) {
                    const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(metadata_pxy_v.at(j));
                    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
                    if (!pfParticlePropertiesMap.empty()) {
                        _topological_score = pfParticlePropertiesMap.at("NuScore");
                    }
                }
            }

            double xyz[3] = {};
            auto vtx = pfp.get<recob::Vertex>();
            if (vtx.size() == 1) {
                vtx.at(0)->XYZ(xyz);
                auto nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);
                _reco_neutrino_vertex_x = nuvtx.X();
                _reco_neutrino_vertex_y = nuvtx.Y();
                _reco_neutrino_vertex_z = nuvtx.Z();
                float reco_nu_vtx_sce[3];
                common::ApplySCECorrectionXYZ(_reco_neutrino_vertex_x, _reco_neutrino_vertex_y, _reco_neutrino_vertex_z, reco_nu_vtx_sce);
                _reco_neutrino_vertex_sce_x = reco_nu_vtx_sce[0];
                _reco_neutrino_vertex_sce_y = reco_nu_vtx_sce[1];
                _reco_neutrino_vertex_sce_z = reco_nu_vtx_sce[2];
            }

            continue;
        }

        _num_pfps++;
        _pfp_slice_indices.push_back(pfpidx++);
        _pfp_pdg_codes.push_back(pfp->PdgCode());
        _pfp_generations.push_back(larpandora.GetGeneration(particleMap, particleMap.at(pfp->Self())));
        uint num_track_daughters = 0;
        uint num_shower_daughters = 0;
        for (size_t daughter : pfp->Daughters()) {
            if (larpandora.IsTrack(particleMap.at(daughter))) {
                num_track_daughters++;
            }
            else {
                num_shower_daughters++;
            }
        }

        _pfp_shower_daughters.push_back(num_shower_daughters);
        _pfp_track_daughters.push_back(num_track_daughters);
        _pfp_num_descendents.push_back(common::GetNDescendents(particleMap.at(pfp->Self()), particleMap));
        
        const auto vertices = pfp.get<recob::Vertex>();
        if (vertices.size() == 1) {
            _pfp_vertex_x.push_back(vertices.at(0)->position().X());
            _pfp_vertex_y.push_back(vertices.at(0)->position().Y());
            _pfp_vertex_z.push_back(vertices.at(0)->position().Z());
        }
        else {
            _pfp_vertex_x.push_back(std::numeric_limits<float>::lowest());
            _pfp_vertex_y.push_back(std::numeric_limits<float>::lowest());
            _pfp_vertex_z.push_back(std::numeric_limits<float>::lowest());
        }
        
        float track_shower_score = common::GetTrackShowerScore(pfp);
        if ((track_shower_score >= 0) && (track_shower_score >= fTrkShrScore)) {
            _num_tracks++;
        }
        else if ((track_shower_score >= 0) && (track_shower_score < fTrkShrScore)) {
            _num_showers++;
        }
        _track_shower_scores.push_back(track_shower_score);

        std::vector<art::Ptr<recob::Hit>> hit_v;
        auto clus_pxy_v = pfp.get<recob::Cluster>();
        _pfp_num_plane_hits_U.push_back(0);
        _pfp_num_plane_hits_V.push_back(0);
        _pfp_num_plane_hits_Y.push_back(0);
        _pfp_num_subclusters_U.push_back(0);
        _pfp_num_subclusters_V.push_back(0);
        _pfp_num_subclusters_Y.push_back(0);
        _pfp_max_subhit_fraction_U.push_back(0);
        _pfp_max_subhit_fraction_V.push_back(0);
        _pfp_max_subhit_fraction_Y.push_back(0);
        for (auto ass_clus : clus_pxy_v) {
            const auto &clus = clus_proxy[ass_clus.key()];
            auto clus_hit_v = clus.get<recob::Hit>();
            auto num_hits = clus_hit_v.size();
            std::vector<art::Ptr<recob::Hit>> cluster_hits_v;
            for (size_t h = 0; h < clus_hit_v.size(); h++) {
                cluster_hits_v.push_back(clus_hit_v[h]);
            }
            int num_subclusters = 0;
            float max_subhit_fraction = 0.;
            std::vector<std::vector<unsigned int>> out_cluster_v;
            if (num_hits) {
                common::cluster(cluster_hits_v, out_cluster_v, 2.0, 1.0);
                for (size_t nc = 0; nc < out_cluster_v.size(); nc++) {
                    auto clus_hit_idx_v = out_cluster_v.at(nc);
                    int num_hits_in_subcluster = clus_hit_idx_v.size();
                    if (num_hits_in_subcluster > 3.) {
                        num_subclusters += 1;
                    }
                    float subhit_fraction = num_hits_in_subcluster / static_cast<float>(num_hits);
                    if (subhit_fraction > max_subhit_fraction) {
                        max_subhit_fraction = subhit_fraction;
                    }
                }
            }
            if (clus->Plane().Plane == 0) {
                _total_hits_U += num_hits;
                _pfp_num_plane_hits_U.back() += num_hits;
                _pfp_num_subclusters_U.back() += num_subclusters;
                _pfp_max_subhit_fraction_U.back() = max_subhit_fraction;
            }
            else if (clus->Plane().Plane == 1) {
                _total_hits_V += num_hits;
                _pfp_num_plane_hits_V.back() += num_hits;
                _pfp_num_subclusters_V.back() += num_subclusters;
                _pfp_max_subhit_fraction_V.back() = max_subhit_fraction;
            }
            else if (clus->Plane().Plane == 2) {
                _total_hits_Y += num_hits;
                _pfp_num_plane_hits_Y.back() += num_hits;
                _pfp_num_subclusters_Y.back() += num_subclusters;
                _pfp_max_subhit_fraction_Y.back() = max_subhit_fraction;
            }
            for (const auto &hit : clus_hit_v) {
                hit_v.push_back(hit);
            }
        }
        _pfp_num_hits.push_back(hit_v.size());

        if (!is_data) {
            if (clus_pxy_v.size() != 0) {
                float purity = 0., completeness = 0., overlay_purity = 0.;
                int ibt = common::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness, overlay_purity);
                if (ibt >= 0) {
                    auto &mcp = btparts_v[ibt];
                    auto PDG = mcp.pdg;
                    _backtracked_energies.push_back(mcp.e);
                    _backtracked_track_ids.push_back(mcp.tids.at(0));
                    _backtracked_pdg_codes.push_back(PDG);
                    _backtracked_purities.push_back(purity);
                    _backtracked_completenesses.push_back(completeness);
                    _backtracked_overlay_purities.push_back(overlay_purity);
                    _backtracked_momentum_x.push_back(mcp.px);
                    _backtracked_momentum_y.push_back(mcp.py);
                    _backtracked_momentum_z.push_back(mcp.pz);
                    _backtracked_start_x.push_back(mcp.start_x);
                    _backtracked_start_y.push_back(mcp.start_y);
                    _backtracked_start_z.push_back(mcp.start_z);
                    _backtracked_start_time.push_back(mcp.start_t);
                    _backtracked_start_wire_U.push_back(common::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 0));
                    _backtracked_start_wire_V.push_back(common::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 1));
                    _backtracked_start_wire_Y.push_back(common::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 2));
                    float reco_st[3] = {mcp.start_x, mcp.start_y, mcp.start_z};
                    if (PDG == 11 || PDG == 22) {
                        reco_st[0] += common::x_offset(mcp.start_t);
                    }
                    else {
                        common::True2RecoMappingXYZ(mcp.start_t, mcp.start_x, mcp.start_y, mcp.start_z, reco_st);
                    }
                    _backtracked_sce_start_x.push_back(reco_st[0]);
                    _backtracked_sce_start_y.push_back(reco_st[1]);
                    _backtracked_sce_start_z.push_back(reco_st[2]);
                    _backtracked_sce_start_wire_U.push_back(common::YZtoPlanecoordinate(reco_st[1], reco_st[2], 0));
                    _backtracked_sce_start_wire_V.push_back(common::YZtoPlanecoordinate(reco_st[1], reco_st[2], 1));
                    _backtracked_sce_start_wire_Y.push_back(common::YZtoPlanecoordinate(reco_st[1], reco_st[2], 2));
                }
                else {
                    _backtracked_energies.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_track_ids.push_back(std::numeric_limits<int>::lowest());
                    _backtracked_pdg_codes.push_back(0);
                    _backtracked_purities.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_completenesses.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_overlay_purities.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_momentum_x.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_momentum_y.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_momentum_z.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_start_x.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_start_y.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_start_z.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_start_time.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_start_wire_U.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_start_wire_V.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_start_wire_Y.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_sce_start_x.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_sce_start_y.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_sce_start_z.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_sce_start_wire_U.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_sce_start_wire_V.push_back(std::numeric_limits<float>::lowest());
                    _backtracked_sce_start_wire_Y.push_back(std::numeric_limits<float>::lowest());
                }
            }
        }
    }

    if (_slice_num_hits > 0) {
        _slice_cluster_fraction = 0.;
        for (auto n : _pfp_num_hits) {
            _slice_cluster_fraction += n;
        }
        _slice_cluster_fraction /= static_cast<float>(_slice_num_hits);
    }

    _num_slices += 1;

    if (is_selected) {
        _selection_pass = 1;
    }
}

void DefaultAnalysis::setBranches(TTree *_tree) {
    _tree->Branch("optical_filter_pe_beam", &_optical_filter_pe_beam, "optical_filter_pe_beam/F");
    _tree->Branch("optical_filter_pe_veto", &_optical_filter_pe_veto, "optical_filter_pe_veto/F");
    _tree->Branch("reco_neutrino_vertex_x", &_reco_neutrino_vertex_x, "reco_neutrino_vertex_x/F");
    _tree->Branch("reco_neutrino_vertex_y", &_reco_neutrino_vertex_y, "reco_neutrino_vertex_y/F");
    _tree->Branch("reco_neutrino_vertex_z", &_reco_neutrino_vertex_z, "reco_neutrino_vertex_z/F");
    _tree->Branch("reco_neutrino_vertex_sce_x", &_reco_neutrino_vertex_sce_x, "reco_neutrino_vertex_sce_x/F");
    _tree->Branch("reco_neutrino_vertex_sce_y", &_reco_neutrino_vertex_sce_y, "reco_neutrino_vertex_sce_y/F");
    _tree->Branch("reco_neutrino_vertex_sce_z", &_reco_neutrino_vertex_sce_z, "reco_neutrino_vertex_sce_z/F");
    _tree->Branch("num_slices", &_num_slices, "num_slices/I");
    _tree->Branch("crt_veto", &_crt_veto, "crt_veto/I");
    _tree->Branch("crt_hit_pe", &_crt_hit_pe, "crt_hit_pe/F");
    _tree->Branch("pfp_slice_indices", "std::vector<int>", &_pfp_slice_indices);
    _tree->Branch("backtracked_pdg_codes", "std::vector<int>", &_backtracked_pdg_codes);
    _tree->Branch("backtracked_energies", "std::vector<float>", &_backtracked_energies);
    _tree->Branch("backtracked_track_ids", "std::vector<int>", &_backtracked_track_ids);
    _tree->Branch("backtracked_purities", "std::vector<float>", &_backtracked_purities);
    _tree->Branch("backtracked_completenesses", "std::vector<float>", &_backtracked_completenesses);
    _tree->Branch("backtracked_overlay_purities", "std::vector<float>", &_backtracked_overlay_purities);
    _tree->Branch("backtracked_momentum_x", "std::vector<float>", &_backtracked_momentum_x);
    _tree->Branch("backtracked_momentum_y", "std::vector<float>", &_backtracked_momentum_y);
    _tree->Branch("backtracked_momentum_z", "std::vector<float>", &_backtracked_momentum_z);
    _tree->Branch("backtracked_start_x", "std::vector<float>", &_backtracked_start_x);
    _tree->Branch("backtracked_start_y", "std::vector<float>", &_backtracked_start_y);
    _tree->Branch("backtracked_start_z", "std::vector<float>", &_backtracked_start_z);
    _tree->Branch("backtracked_start_time", "std::vector<float>", &_backtracked_start_time);
    _tree->Branch("backtracked_start_wire_U", "std::vector<float>", &_backtracked_start_wire_U);
    _tree->Branch("backtracked_start_wire_V", "std::vector<float>", &_backtracked_start_wire_V);
    _tree->Branch("backtracked_start_wire_Y", "std::vector<float>", &_backtracked_start_wire_Y);
    _tree->Branch("backtracked_sce_start_x", "std::vector<float>", &_backtracked_sce_start_x);
    _tree->Branch("backtracked_sce_start_y", "std::vector<float>", &_backtracked_sce_start_y);
    _tree->Branch("backtracked_sce_start_z", "std::vector<float>", &_backtracked_sce_start_z);
    _tree->Branch("backtracked_sce_start_wire_U", "std::vector<float>", &_backtracked_sce_start_wire_U);
    _tree->Branch("backtracked_sce_start_wire_V", "std::vector<float>", &_backtracked_sce_start_wire_V);
    _tree->Branch("backtracked_sce_start_wire_Y", "std::vector<float>", &_backtracked_sce_start_wire_Y);
    _tree->Branch("selection_pass", &_selection_pass, "selection_pass/I");
    _tree->Branch("software_trigger", &_software_trigger, "software_trigger/I");
    if (fMakeNuMINtuple) {
        _tree->Branch("software_trigger_pre", &_software_trigger_pre, "software_trigger_pre/I");
        _tree->Branch("software_trigger_post", &_software_trigger_post, "software_trigger_post/I");
        _tree->Branch("software_trigger_pre_ext", &_software_trigger_pre_ext, "software_trigger_pre_ext/I");
        _tree->Branch("software_trigger_post_ext", &_software_trigger_post_ext, "software_trigger_post_ext/I");
    }
    _tree->Branch("event_total_hits", &_event_total_hits, "event_total_hits/I");
    _tree->Branch("slice_pdg", &_slice_pdg, "slice_pdg/I");
    _tree->Branch("slice_num_hits", &_slice_num_hits, "slice_num_hits/I");
    _tree->Branch("num_pfps", &_num_pfps, "num_pfps/I");
    _tree->Branch("num_tracks", &_num_tracks, "num_tracks/I");
    _tree->Branch("num_showers", &_num_showers, "num_showers/I");
    _tree->Branch("pfp_generations", "std::vector<uint>", &_pfp_generations);
    _tree->Branch("pfp_track_daughters", "std::vector<uint>", &_pfp_track_daughters);
    _tree->Branch("pfp_shower_daughters", "std::vector<uint>", &_pfp_shower_daughters);
    _tree->Branch("pfp_num_descendents", "std::vector<uint>", &_pfp_num_descendents);
    _tree->Branch("pfp_vertex_x", "std::vector<float>", &_pfp_vertex_x);
    _tree->Branch("pfp_vertex_y", "std::vector<float>", &_pfp_vertex_y);
    _tree->Branch("pfp_vertex_z", "std::vector<float>", &_pfp_vertex_z);
    _tree->Branch("track_shower_scores", "std::vector<float>", &_track_shower_scores);
    _tree->Branch("pfp_pdg_codes", "std::vector<int>", &_pfp_pdg_codes);
    _tree->Branch("pfp_num_hits", "std::vector<int>", &_pfp_num_hits);
    _tree->Branch("pfp_num_plane_hits_U", "std::vector<int>", &_pfp_num_plane_hits_U);
    _tree->Branch("pfp_num_plane_hits_V", "std::vector<int>", &_pfp_num_plane_hits_V);
    _tree->Branch("pfp_num_plane_hits_Y", "std::vector<int>", &_pfp_num_plane_hits_Y);
    _tree->Branch("pfp_num_subclusters_U", "std::vector<int>", &_pfp_num_subclusters_U);
    _tree->Branch("pfp_num_subclusters_V", "std::vector<int>", &_pfp_num_subclusters_V);
    _tree->Branch("pfp_num_subclusters_Y", "std::vector<int>", &_pfp_num_subclusters_Y);
    _tree->Branch("pfp_max_subhit_fraction_U", "std::vector<float>", &_pfp_max_subhit_fraction_U);
    _tree->Branch("pfp_max_subhit_fraction_V", "std::vector<float>", &_pfp_max_subhit_fraction_V);
    _tree->Branch("pfp_max_subhit_fraction_Y", "std::vector<float>", &_pfp_max_subhit_fraction_Y);
    _tree->Branch("total_hits_U", &_total_hits_U, "total_hits_U/i");
    _tree->Branch("total_hits_V", &_total_hits_V, "total_hits_V/i");
    _tree->Branch("total_hits_Y", &_total_hits_Y, "total_hits_Y/i");
    _tree->Branch("slice_id", &_slice_id, "slice_id/i");
    _tree->Branch("slice_topological_scores", "std::vector<float>", &_slice_topological_scores);
    _tree->Branch("topological_score", &_topological_score, "topological_score/F");
    _tree->Branch("slice_cluster_fraction", &_slice_cluster_fraction, "slice_cluster_fraction/F");
}

void DefaultAnalysis::resetTTree(TTree *_tree) {
    if (fMakeNuMINtuple) {
        _software_trigger = 1;
        _software_trigger_pre = 1;
        _software_trigger_post = 1;
        _software_trigger_pre_ext = 1;
        _software_trigger_post_ext = 1;
    }
    else {
        _software_trigger = 0;
    }

    _selection_pass = 0;
    _optical_filter_pe_beam = 0.;
    _optical_filter_pe_veto = 0.;
    _reco_neutrino_vertex_x = std::numeric_limits<float>::lowest();
    _reco_neutrino_vertex_y = std::numeric_limits<float>::lowest();
    _reco_neutrino_vertex_z = std::numeric_limits<float>::lowest();
    _reco_neutrino_vertex_sce_x = std::numeric_limits<float>::lowest();
    _reco_neutrino_vertex_sce_y = std::numeric_limits<float>::lowest();
    _reco_neutrino_vertex_sce_z = std::numeric_limits<float>::lowest();
    _num_slices = 0;
    _crt_veto = 0;
    _crt_hit_pe = 0;
    _num_pfps = 0;
    _num_tracks = 0;
    _num_showers = 0;

    _track_shower_scores.clear();
    _backtracked_energies.clear();
    _backtracked_track_ids.clear();
    _backtracked_pdg_codes.clear();
    _backtracked_purities.clear();
    _backtracked_completenesses.clear();
    _backtracked_overlay_purities.clear();
    _backtracked_momentum_x.clear();
    _backtracked_momentum_y.clear();
    _backtracked_momentum_z.clear();
    _backtracked_start_x.clear();
    _backtracked_start_y.clear();
    _backtracked_start_z.clear();
    _backtracked_start_time.clear();
    _backtracked_start_wire_U.clear();
    _backtracked_start_wire_V.clear();
    _backtracked_start_wire_Y.clear();
    _backtracked_sce_start_x.clear();
    _backtracked_sce_start_y.clear();
    _backtracked_sce_start_z.clear();
    _backtracked_sce_start_wire_U.clear();
    _backtracked_sce_start_wire_V.clear();
    _backtracked_sce_start_wire_Y.clear();

    _event_total_hits = std::numeric_limits<int>::lowest();
    _slice_pdg = std::numeric_limits<int>::lowest();
    _slice_id = std::numeric_limits<int>::lowest();
    _topological_score = std::numeric_limits<float>::lowest();
    _slice_topological_scores.clear();
    _slice_num_hits = std::numeric_limits<int>::lowest();
    _pfp_pdg_codes.clear();
    _pfp_num_hits.clear();
    _pfp_num_plane_hits_U.clear();
    _pfp_num_plane_hits_V.clear();
    _pfp_num_plane_hits_Y.clear();
    _pfp_num_subclusters_U.clear();
    _pfp_num_subclusters_V.clear();
    _pfp_num_subclusters_Y.clear();
    _pfp_max_subhit_fraction_U.clear();
    _pfp_max_subhit_fraction_V.clear();
    _pfp_max_subhit_fraction_Y.clear();
    _pfp_generations.clear();
    _pfp_shower_daughters.clear();
    _pfp_track_daughters.clear();
    _pfp_num_descendents.clear();
    _pfp_vertex_x.clear();
    _pfp_vertex_y.clear();
    _pfp_vertex_z.clear();
    _slice_cluster_fraction = std::numeric_limits<float>::lowest();
    _total_hits_U = 0;
    _total_hits_V = 0;
    _total_hits_Y = 0;
}

DEFINE_ART_CLASS_TOOL(DefaultAnalysis)


}

#endif
