#ifndef DEFAULT_ANALYSIS_CXX
#define DEFAULT_ANALYSIS_CXX

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"

#include "AnalysisToolBase.h"
#include "Helpers/GeometryUtils.h"
#include "Helpers/HitProximityClustering.h"
#include "Helpers/PFParticleDescendants.h"
#include "Helpers/PandoraUtilities.h"
#include "Helpers/SpaceChargeCorrections.h"
#include "Helpers/TrackShowerScore.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>

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
    art::InputTag fCLSproducer;
    art::InputTag fMCTproducer;
    art::InputTag fMCPproducer;
    art::InputTag fMCFproducer;
    art::InputTag fHITproducer;
    art::InputTag fSLCproducer;
    art::InputTag fSPproducer;

    float fTrkShrScore;
    std::string NuMIOpFilterProd;
    std::string NuMISWTrigProd;

    bool fMakeNuMINtuple;
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

    std::vector<int> _pfp_slice_indices;
    int _selection_pass;
    int _event_total_hits;
    int _slice_pdg;
    int _slice_num_hits;
    int _slice_num_hits_u;
    int _slice_num_hits_v;
    int _slice_num_hits_w;
    int _slice_num_clusters;
    int _slice_num_clusters_u;
    int _slice_num_clusters_v;
    int _slice_num_clusters_w;
    int _slice_id;
    float _topological_score;

    std::vector<int> _pfp_pdg_codes;
    std::vector<int> _pfp_num_hits;
    std::vector<int> _pfp_num_plane_hits_u;
    std::vector<int> _pfp_num_plane_hits_v;
    std::vector<int> _pfp_num_plane_hits_w;
    int _slice_num_spacepoints;
    std::vector<float> _slice_spacepoint_x;
    std::vector<float> _slice_spacepoint_y;
    std::vector<float> _slice_spacepoint_z;
    std::vector<float> _slice_spacepoint_charge;
    std::vector<float> _slice_spacepoint_extent_x;
    std::vector<float> _slice_spacepoint_extent_y;
    std::vector<float> _slice_spacepoint_extent_z;
    std::vector<float> _slice_spacepoint_extent_3d;
    std::vector<float> _slice_hit_charge_u;
    std::vector<float> _slice_hit_charge_v;
    std::vector<float> _slice_hit_charge_w;
    std::vector<int> _slice_hit_wire_u;
    std::vector<int> _slice_hit_wire_v;
    std::vector<int> _slice_hit_wire_w;
    std::vector<float> _slice_hit_peak_time_u;
    std::vector<float> _slice_hit_peak_time_v;
    std::vector<float> _slice_hit_peak_time_w;
    std::vector<int> _slice_hit_start_tick_u;
    std::vector<int> _slice_hit_start_tick_v;
    std::vector<int> _slice_hit_start_tick_w;
    std::vector<int> _slice_hit_end_tick_u;
    std::vector<int> _slice_hit_end_tick_v;
    std::vector<int> _slice_hit_end_tick_w;
    std::vector<int> _slice_hit_tick_extent_u;
    std::vector<int> _slice_hit_tick_extent_v;
    std::vector<int> _slice_hit_tick_extent_w;
    std::vector<int> _pfp_num_subclusters_u;
    std::vector<int> _pfp_num_subclusters_v;
    std::vector<int> _pfp_num_subclusters_w;
    std::vector<float> _pfp_max_subhit_fraction_u;
    std::vector<float> _pfp_max_subhit_fraction_v;
    std::vector<float> _pfp_max_subhit_fraction_w;

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
    unsigned int _total_hits_u;
    unsigned int _total_hits_v;
    unsigned int _total_hits_w;
};

DefaultAnalysis::DefaultAnalysis(const fhicl::ParameterSet &p) {
    fPFPproducer = p.get<art::InputTag>("PFPproducer", "pandora");
    fCLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
    fMCTproducer = p.get<art::InputTag>("MCTproducer", "generator");
    fMCPproducer = p.get<art::InputTag>("MCPproducer", "largeant");
    fMCFproducer = p.get<art::InputTag>("MCFproducer", "generator");
    fHITproducer = p.get<art::InputTag>("HITproducer", "gaushit");
    fSLCproducer = p.get<art::InputTag>("SLCproducer", "pandora");
    fSPproducer = p.get<art::InputTag>("SPproducer", "pandora");

    fTrkShrScore = p.get<float>("TrkShrScore", 0.5);
    
    fMakeNuMINtuple = p.get<bool>("makeNuMINtuple", true);
    NuMIOpFilterProd = p.get<std::string>("NuMIOpFiltProcName", "DataStage1Optical");
    NuMISWTrigProd = p.get<std::string>("NuMISWTriggerProcName", "DataOverlayOpticalNuMI");
}

void DefaultAnalysis::configure(fhicl::ParameterSet const &p) {
}

void DefaultAnalysis::analyseEvent(const art::Event &event, bool is_data) {
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
            auto const& trigger = *triggerHandle;
            for (int j = 0; j != triggerHandle->getNumberOfAlgorithms(); j++) {
                const bool passed = trigger.passedAlgo(triggerName[j]);
                if (triggerName[j] == "EXT_NUMIwin_FEMBeamTriggerAlgo") {
                    _software_trigger_pre_ext = passed;
                }
                else if (triggerName[j] == "EXT_NUMIwin_2018May_FEMBeamTriggerAlgo") {
                    _software_trigger_post_ext = passed;
                }
                else if (triggerName[j] == "NUMI_FEMBeamTriggerAlgo") {
                    _software_trigger_pre = passed;
                }
                else if (triggerName[j] == "NUMI_2018May_FEMBeamTriggerAlgo") {
                    _software_trigger_post = passed;
                }
                else {
                    continue;
                }
                std::cout << triggerName[j] << ": ";
                std::cout << passed << std::endl;
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

    art::ValidHandle<std::vector<recob::Hit>> inputHits = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    _event_total_hits = inputHits->size();
}

void DefaultAnalysis::analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool is_data, bool is_selected) {
    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(event, fCLSproducer,
        proxy::withAssociated<recob::Hit>(fCLSproducer));
    art::ValidHandle<std::vector<recob::Slice>> inputSlice = event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, event, fSLCproducer));
    art::FindManyP<recob::SpacePoint> hit_to_sp(hit_handle, event, fSPproducer);

    lar_pandora::LArPandoraHelper larpandora;
    lar_pandora::PFParticleVector pfparticles;
    lar_pandora::PFParticleMap particleMap;
    larpandora.CollectPFParticles(event, "pandora", pfparticles);
    larpandora.BuildPFParticleMap(pfparticles, particleMap);

    size_t pfpidx = 0;
    _num_pfps = 0;
    _slice_num_hits_u = 0;
    _slice_num_hits_v = 0;
    _slice_num_hits_w = 0;
    _slice_num_spacepoints = 0;
    _slice_spacepoint_x.clear();
    _slice_spacepoint_y.clear();
    _slice_spacepoint_z.clear();
    _slice_spacepoint_charge.clear();
    _slice_spacepoint_extent_x.clear();
    _slice_spacepoint_extent_y.clear();
    _slice_spacepoint_extent_z.clear();
    _slice_spacepoint_extent_3d.clear();
    _slice_hit_charge_u.clear();
    _slice_hit_charge_v.clear();
    _slice_hit_charge_w.clear();
    _slice_hit_wire_u.clear();
    _slice_hit_wire_v.clear();
    _slice_hit_wire_w.clear();
    _slice_hit_peak_time_u.clear();
    _slice_hit_peak_time_v.clear();
    _slice_hit_peak_time_w.clear();
    _slice_hit_start_tick_u.clear();
    _slice_hit_start_tick_v.clear();
    _slice_hit_start_tick_w.clear();
    _slice_hit_end_tick_u.clear();
    _slice_hit_end_tick_v.clear();
    _slice_hit_end_tick_w.clear();
    _slice_hit_tick_extent_u.clear();
    _slice_hit_tick_extent_v.clear();
    _slice_hit_tick_extent_w.clear();
    _slice_num_clusters = 0;
    _slice_num_clusters_u = 0;
    _slice_num_clusters_v = 0;
    _slice_num_clusters_w = 0;
    for (auto pfp : slice_pfp_vec) {
        if (pfp->IsPrimary()) {
            _slice_pdg = pfp->PdgCode();
            auto slice_pxy_v = pfp.get<recob::Slice>();
            if (slice_pxy_v.size() != 1) 
                return;
            
            auto slicehits = assocSliceHit->at(slice_pxy_v[0].key());
            _slice_num_hits = slicehits.size();
            std::map<art::Ptr<recob::SpacePoint>, double> slice_spacepoint_charge_map;
            std::set<std::size_t> seen_spacepoint_keys;
            auto add_slice_hit = [&](const art::Ptr<recob::Hit> &hit) {
                std::vector<float> *charges = nullptr;
                std::vector<int> *wires = nullptr;
                std::vector<float> *peak_times = nullptr;
                std::vector<int> *start_ticks = nullptr;
                std::vector<int> *end_ticks = nullptr;
                std::vector<int> *tick_extents = nullptr;
                int *plane_hit_count = nullptr;

                if (hit->WireID().Plane == 0) {
                    charges = &_slice_hit_charge_u;
                    wires = &_slice_hit_wire_u;
                    peak_times = &_slice_hit_peak_time_u;
                    start_ticks = &_slice_hit_start_tick_u;
                    end_ticks = &_slice_hit_end_tick_u;
                    tick_extents = &_slice_hit_tick_extent_u;
                    plane_hit_count = &_slice_num_hits_u;
                }
                else if (hit->WireID().Plane == 1) {
                    charges = &_slice_hit_charge_v;
                    wires = &_slice_hit_wire_v;
                    peak_times = &_slice_hit_peak_time_v;
                    start_ticks = &_slice_hit_start_tick_v;
                    end_ticks = &_slice_hit_end_tick_v;
                    tick_extents = &_slice_hit_tick_extent_v;
                    plane_hit_count = &_slice_num_hits_v;
                }
                else if (hit->WireID().Plane == 2) {
                    charges = &_slice_hit_charge_w;
                    wires = &_slice_hit_wire_w;
                    peak_times = &_slice_hit_peak_time_w;
                    start_ticks = &_slice_hit_start_tick_w;
                    end_ticks = &_slice_hit_end_tick_w;
                    tick_extents = &_slice_hit_tick_extent_w;
                    plane_hit_count = &_slice_num_hits_w;
                }
                else {
                    return;
                }

                *plane_hit_count += 1;
                charges->push_back(hit->Integral());
                wires->push_back(static_cast<int>(hit->WireID().Wire));
                peak_times->push_back(static_cast<float>(hit->PeakTime()));
                start_ticks->push_back(static_cast<int>(hit->StartTick()));
                end_ticks->push_back(static_cast<int>(hit->EndTick()));
                tick_extents->push_back(static_cast<int>(hit->EndTick() - hit->StartTick()));

                double const hit_charge = std::max(0.f, hit->Integral());
                if (!(hit_charge > 0.0))
                    return;

                auto const &spacepoints = hit_to_sp.at(hit.key());
                if (spacepoints.empty())
                    return;

                double const charge_per_spacepoint = hit_charge / static_cast<double>(spacepoints.size());
                for (auto const &spacepoint : spacepoints) {
                    if (!spacepoint)
                        continue;
                    slice_spacepoint_charge_map[spacepoint] += charge_per_spacepoint;
                }
            };
            for (auto const &hit : slicehits) {
                add_slice_hit(hit);
            }
            for (auto const &pfp_in_slice : slice_pfp_vec) {
                auto const &spacepoints = pfp_in_slice.get<recob::SpacePoint>();
                for (auto const &spacepoint : spacepoints) {
                    if (!spacepoint)
                        continue;
                    if (!seen_spacepoint_keys.insert(static_cast<std::size_t>(spacepoint.key())).second)
                        continue;

                    auto const *xyz = spacepoint->XYZ();
                    auto const *err = spacepoint->ErrXYZ();
                    float const extent_x = static_cast<float>(std::sqrt(std::max(0., static_cast<double>(err[0]))));
                    float const extent_y = static_cast<float>(std::sqrt(std::max(0., static_cast<double>(err[2]))));
                    float const extent_z = static_cast<float>(std::sqrt(std::max(0., static_cast<double>(err[5]))));

                    _slice_spacepoint_x.push_back(static_cast<float>(xyz[0]));
                    _slice_spacepoint_y.push_back(static_cast<float>(xyz[1]));
                    _slice_spacepoint_z.push_back(static_cast<float>(xyz[2]));
                    _slice_spacepoint_charge.push_back(static_cast<float>(slice_spacepoint_charge_map[spacepoint]));
                    _slice_spacepoint_extent_x.push_back(extent_x);
                    _slice_spacepoint_extent_y.push_back(extent_y);
                    _slice_spacepoint_extent_z.push_back(extent_z);
                    _slice_spacepoint_extent_3d.push_back(std::sqrt(extent_x * extent_x + extent_y * extent_y + extent_z * extent_z));
                }
            }
            _slice_num_spacepoints = static_cast<int>(_slice_spacepoint_x.size());
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
            _pfp_vertex_x.push_back(std::numeric_limits<float>::quiet_NaN());
            _pfp_vertex_y.push_back(std::numeric_limits<float>::quiet_NaN());
            _pfp_vertex_z.push_back(std::numeric_limits<float>::quiet_NaN());
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
        _pfp_num_plane_hits_u.push_back(0);
        _pfp_num_plane_hits_v.push_back(0);
        _pfp_num_plane_hits_w.push_back(0);
        _pfp_num_subclusters_u.push_back(0);
        _pfp_num_subclusters_v.push_back(0);
        _pfp_num_subclusters_w.push_back(0);
        _pfp_max_subhit_fraction_u.push_back(0);
        _pfp_max_subhit_fraction_v.push_back(0);
        _pfp_max_subhit_fraction_w.push_back(0);
        for (auto ass_clus : clus_pxy_v) {
            const auto &clus = clus_proxy[ass_clus.key()];
            _slice_num_clusters += 1;
            if (clus->Plane().Plane == 0) {
                _slice_num_clusters_u += 1;
            }
            else if (clus->Plane().Plane == 1) {
                _slice_num_clusters_v += 1;
            }
            else if (clus->Plane().Plane == 2) {
                _slice_num_clusters_w += 1;
            }
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
                _total_hits_u += num_hits;
                _pfp_num_plane_hits_u.back() += num_hits;
                _pfp_num_subclusters_u.back() += num_subclusters;
                _pfp_max_subhit_fraction_u.back() = max_subhit_fraction;
            }
            else if (clus->Plane().Plane == 1) {
                _total_hits_v += num_hits;
                _pfp_num_plane_hits_v.back() += num_hits;
                _pfp_num_subclusters_v.back() += num_subclusters;
                _pfp_max_subhit_fraction_v.back() = max_subhit_fraction;
            }
            else if (clus->Plane().Plane == 2) {
                _total_hits_w += num_hits;
                _pfp_num_plane_hits_w.back() += num_hits;
                _pfp_num_subclusters_w.back() += num_subclusters;
                _pfp_max_subhit_fraction_w.back() = max_subhit_fraction;
            }
            for (const auto &hit : clus_hit_v) {
                hit_v.push_back(hit);
            }
        }
        _pfp_num_hits.push_back(hit_v.size());

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
    _tree->Branch("pfp_slice_indices", "std::vector<int>", &_pfp_slice_indices);
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
    _tree->Branch("slice_num_hits_u", &_slice_num_hits_u, "slice_num_hits_u/I");
    _tree->Branch("slice_num_hits_v", &_slice_num_hits_v, "slice_num_hits_v/I");
    _tree->Branch("slice_num_hits_w", &_slice_num_hits_w, "slice_num_hits_w/I");
    _tree->Branch("slice_num_spacepoints", &_slice_num_spacepoints, "slice_num_spacepoints/I");
    _tree->Branch("slice_num_clusters", &_slice_num_clusters, "slice_num_clusters/I");
    _tree->Branch("slice_num_clusters_u", &_slice_num_clusters_u, "slice_num_clusters_u/I");
    _tree->Branch("slice_num_clusters_v", &_slice_num_clusters_v, "slice_num_clusters_v/I");
    _tree->Branch("slice_num_clusters_w", &_slice_num_clusters_w, "slice_num_clusters_w/I");
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
    _tree->Branch("pfp_num_plane_hits_u", "std::vector<int>", &_pfp_num_plane_hits_u);
    _tree->Branch("pfp_num_plane_hits_v", "std::vector<int>", &_pfp_num_plane_hits_v);
    _tree->Branch("pfp_num_plane_hits_w", "std::vector<int>", &_pfp_num_plane_hits_w);
    _tree->Branch("slice_spacepoint_x", "std::vector<float>", &_slice_spacepoint_x);
    _tree->Branch("slice_spacepoint_y", "std::vector<float>", &_slice_spacepoint_y);
    _tree->Branch("slice_spacepoint_z", "std::vector<float>", &_slice_spacepoint_z);
    _tree->Branch("slice_spacepoint_charge", "std::vector<float>", &_slice_spacepoint_charge);
    _tree->Branch("slice_spacepoint_extent_x", "std::vector<float>", &_slice_spacepoint_extent_x);
    _tree->Branch("slice_spacepoint_extent_y", "std::vector<float>", &_slice_spacepoint_extent_y);
    _tree->Branch("slice_spacepoint_extent_z", "std::vector<float>", &_slice_spacepoint_extent_z);
    _tree->Branch("slice_spacepoint_extent_3d", "std::vector<float>", &_slice_spacepoint_extent_3d);
    _tree->Branch("slice_hit_charge_u", "std::vector<float>", &_slice_hit_charge_u);
    _tree->Branch("slice_hit_charge_v", "std::vector<float>", &_slice_hit_charge_v);
    _tree->Branch("slice_hit_charge_w", "std::vector<float>", &_slice_hit_charge_w);
    _tree->Branch("slice_hit_wire_u", "std::vector<int>", &_slice_hit_wire_u);
    _tree->Branch("slice_hit_wire_v", "std::vector<int>", &_slice_hit_wire_v);
    _tree->Branch("slice_hit_wire_w", "std::vector<int>", &_slice_hit_wire_w);
    _tree->Branch("slice_hit_peak_time_u", "std::vector<float>", &_slice_hit_peak_time_u);
    _tree->Branch("slice_hit_peak_time_v", "std::vector<float>", &_slice_hit_peak_time_v);
    _tree->Branch("slice_hit_peak_time_w", "std::vector<float>", &_slice_hit_peak_time_w);
    _tree->Branch("slice_hit_start_tick_u", "std::vector<int>", &_slice_hit_start_tick_u);
    _tree->Branch("slice_hit_start_tick_v", "std::vector<int>", &_slice_hit_start_tick_v);
    _tree->Branch("slice_hit_start_tick_w", "std::vector<int>", &_slice_hit_start_tick_w);
    _tree->Branch("slice_hit_end_tick_u", "std::vector<int>", &_slice_hit_end_tick_u);
    _tree->Branch("slice_hit_end_tick_v", "std::vector<int>", &_slice_hit_end_tick_v);
    _tree->Branch("slice_hit_end_tick_w", "std::vector<int>", &_slice_hit_end_tick_w);
    _tree->Branch("slice_hit_tick_extent_u", "std::vector<int>", &_slice_hit_tick_extent_u);
    _tree->Branch("slice_hit_tick_extent_v", "std::vector<int>", &_slice_hit_tick_extent_v);
    _tree->Branch("slice_hit_tick_extent_w", "std::vector<int>", &_slice_hit_tick_extent_w);
    _tree->Branch("pfp_num_subclusters_u", "std::vector<int>", &_pfp_num_subclusters_u);
    _tree->Branch("pfp_num_subclusters_v", "std::vector<int>", &_pfp_num_subclusters_v);
    _tree->Branch("pfp_num_subclusters_w", "std::vector<int>", &_pfp_num_subclusters_w);
    _tree->Branch("pfp_max_subhit_fraction_u", "std::vector<float>", &_pfp_max_subhit_fraction_u);
    _tree->Branch("pfp_max_subhit_fraction_v", "std::vector<float>", &_pfp_max_subhit_fraction_v);
    _tree->Branch("pfp_max_subhit_fraction_w", "std::vector<float>", &_pfp_max_subhit_fraction_w);
    _tree->Branch("total_hits_u", &_total_hits_u, "total_hits_u/i");
    _tree->Branch("total_hits_v", &_total_hits_v, "total_hits_v/i");
    _tree->Branch("total_hits_w", &_total_hits_w, "total_hits_w/i");
    _tree->Branch("slice_id", &_slice_id, "slice_id/i");
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
    _reco_neutrino_vertex_x = std::numeric_limits<float>::quiet_NaN();
    _reco_neutrino_vertex_y = std::numeric_limits<float>::quiet_NaN();
    _reco_neutrino_vertex_z = std::numeric_limits<float>::quiet_NaN();
    _reco_neutrino_vertex_sce_x = std::numeric_limits<float>::quiet_NaN();
    _reco_neutrino_vertex_sce_y = std::numeric_limits<float>::quiet_NaN();
    _reco_neutrino_vertex_sce_z = std::numeric_limits<float>::quiet_NaN();
    _num_slices = 0;
    _num_pfps = 0;
    _num_tracks = 0;
    _num_showers = 0;

    _track_shower_scores.clear();

    _event_total_hits = -1;
    _slice_pdg = 0;
    _slice_id = -1;
    _topological_score = std::numeric_limits<float>::quiet_NaN();
    _slice_num_hits = -1;
    _slice_num_hits_u = 0;
    _slice_num_hits_v = 0;
    _slice_num_hits_w = 0;
    _slice_num_spacepoints = 0;
    _slice_num_clusters = 0;
    _slice_num_clusters_u = 0;
    _slice_num_clusters_v = 0;
    _slice_num_clusters_w = 0;
    _pfp_slice_indices.clear();
    _pfp_pdg_codes.clear();
    _pfp_num_hits.clear();
    _pfp_num_plane_hits_u.clear();
    _pfp_num_plane_hits_v.clear();
    _pfp_num_plane_hits_w.clear();
    _slice_spacepoint_x.clear();
    _slice_spacepoint_y.clear();
    _slice_spacepoint_z.clear();
    _slice_spacepoint_charge.clear();
    _slice_spacepoint_extent_x.clear();
    _slice_spacepoint_extent_y.clear();
    _slice_spacepoint_extent_z.clear();
    _slice_spacepoint_extent_3d.clear();
    _slice_hit_charge_u.clear();
    _slice_hit_charge_v.clear();
    _slice_hit_charge_w.clear();
    _slice_hit_wire_u.clear();
    _slice_hit_wire_v.clear();
    _slice_hit_wire_w.clear();
    _slice_hit_peak_time_u.clear();
    _slice_hit_peak_time_v.clear();
    _slice_hit_peak_time_w.clear();
    _slice_hit_start_tick_u.clear();
    _slice_hit_start_tick_v.clear();
    _slice_hit_start_tick_w.clear();
    _slice_hit_end_tick_u.clear();
    _slice_hit_end_tick_v.clear();
    _slice_hit_end_tick_w.clear();
    _slice_hit_tick_extent_u.clear();
    _slice_hit_tick_extent_v.clear();
    _slice_hit_tick_extent_w.clear();
    _pfp_num_subclusters_u.clear();
    _pfp_num_subclusters_v.clear();
    _pfp_num_subclusters_w.clear();
    _pfp_max_subhit_fraction_u.clear();
    _pfp_max_subhit_fraction_v.clear();
    _pfp_max_subhit_fraction_w.clear();
    _pfp_generations.clear();
    _pfp_shower_daughters.clear();
    _pfp_track_daughters.clear();
    _pfp_num_descendents.clear();
    _pfp_vertex_x.clear();
    _pfp_vertex_y.clear();
    _pfp_vertex_z.clear();
    _slice_cluster_fraction = std::numeric_limits<float>::quiet_NaN();
    _total_hits_u = 0;
    _total_hits_v = 0;
    _total_hits_w = 0;
}

DEFINE_ART_CLASS_TOOL(DefaultAnalysis)

}

#endif
