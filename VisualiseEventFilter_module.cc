#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/MCBase/MCShower.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "CommonFunctions/Types.h"
#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"

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

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"

#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <tuple>

class VisualiseEventFilter : public art::EDFilter
{
public:
    explicit VisualiseEventFilter(fhicl::ParameterSet const &pset);

    VisualiseEventFilter(VisualiseEventFilter const &) = delete;
    VisualiseEventFilter(VisualiseEventFilter &&) = delete;
    VisualiseEventFilter &operator=(VisualiseEventFilter const &) = delete;
    VisualiseEventFilter &operator=(VisualiseEventFilter &&) = delete;

    bool filter(art::Event &e) override;

private:
    art::InputTag _HitProducer, _MCTproducer, _BacktrackTag;

    std::vector<std::tuple<int, int, int>> _target_events;
    std::string _mode;

    void visualiseEvent(const art::Event &e, const std::string &filename);
    void get_limits(const std::vector<float>& wire_coord_vec, const std::vector<float>& drift_coord_vec, float& global_wire_min, float& global_wire_max, float& global_drift_min, float& global_drift_max) const;
};

VisualiseEventFilter::VisualiseEventFilter(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
    , _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "largeant")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _mode{pset.get<std::string>("Mode", "target")}
{
    bool has_target_events = pset.has_key("TargetEvents");
    bool has_target_events_file = pset.has_key("TargetEventsFile");

    if (has_target_events && has_target_events_file) 
        throw cet::exception("VisualiseEventFilter") << "-- Cannot specify both 'TargetEvents' and 'TargetEventsFile'. Please choose one.";

    if (has_target_events) {
        auto entry_coll = pset.get<std::vector<std::vector<int>>>("TargetEvents");
        for (const auto &entry : entry_coll) {
            if (entry.size() == 3) 
                _target_events.emplace_back(entry[0], entry[1], entry[2]);
            else 
                throw cet::exception("VisualiseEventFilter") << "-- Invalid TargetEvents format. Each entry must contain exactly three integers (Run, SubRun, Event).";
        }
    } else if (has_target_events_file) {
        std::string file_path = pset.get<std::string>("TargetEventsFile");
        std::ifstream file(file_path);
        if (!file.is_open()) 
            throw cet::exception("VisualiseEventFilter") << "-- Unable to open TargetEventsFile: " << file_path;

        int run, subrun, event;
        while (file >> run >> subrun >> event) 
            _target_events.emplace_back(run, subrun, event);
    } else 
        throw cet::exception("VisualiseEventFilter") << "-- Either 'TargetEvents' or 'TargetEventsFile' must be specified.";
}

bool VisualiseEventFilter::filter(art::Event &e)
{
    auto current_event = std::make_tuple(e.run(), e.subRun(), e.event());
    if (_mode == "target" && std::find(_target_events.begin(), _target_events.end(), current_event) == _target_events.end()) 
        return false; 

    mf::LogInfo("VisualiseEventFilter") << "-- Visualising event: Run " << e.run() 
                                           << ", SubRun " << e.subRun() 
                                           << ", Event " << e.event();

    std::string filename = "event_" + std::to_string(e.run()) + "_" +
                           std::to_string(e.subRun()) + "_" +
                           std::to_string(e.event()) + ".png";

    this->visualiseEvent(e, filename);

    return true;
}

void VisualiseEventFilter::visualiseEvent(const art::Event &e, const std::string &filename)
{
    art::Handle<std::vector<simb::MCParticle>> mc_particle_handle; 
    std::vector<art::Ptr<simb::MCParticle>> mc_particle_vector;
    lar_pandora::MCParticleMap mc_particle_map;

    if (!e.getByLabel(_MCTproducer, mc_particle_handle))
        throw cet::exception("VisualiseEventFilter") << "failed to find any mc particles in event" << std::endl;
    art::fill_ptr_vector(mc_particle_vector, mc_particle_handle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mc_particle_vector, mc_particle_map);

    art::Handle<std::vector<recob::Hit>> hit_handle;
    std::vector<art::Ptr<recob::Hit>> hit_vector;
    
    if (!e.getByLabel(_HitProducer, hit_handle))
        throw cet::exception("VisualiseEventFilter") << "failed to find any hits in event" << std::endl;
    art::fill_ptr_vector(hit_vector, hit_handle);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hit_handle, e, _BacktrackTag);

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

    std::vector<float> true_hits_u_wire;
    std::vector<float> true_hits_u_drift;
    std::vector<float> true_hits_u_owner;
    std::vector<float> true_hits_v_wire;
    std::vector<float> true_hits_v_drift;
    std::vector<float> true_hits_v_owner;
    std::vector<float> true_hits_w_wire;
    std::vector<float> true_hits_w_drift;
    std::vector<float> true_hits_w_owner;

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
            true_hits_u_wire.push_back(pandora_pos.Z());
            true_hits_u_drift.push_back(pandora_pos.X());
            true_hits_u_owner.push_back(owner_pdg_code);
        }
        else if (pandora_view == common::TPC_VIEW_V) {
            true_hits_v_wire.push_back(pandora_pos.Z());
            true_hits_v_drift.push_back(pandora_pos.X());
            true_hits_v_owner.push_back(owner_pdg_code);
        }
        else if (pandora_view == common::TPC_VIEW_W) {
            true_hits_w_wire.push_back(pandora_pos.Z());
            true_hits_w_drift.push_back(pandora_pos.X());
            true_hits_w_owner.push_back(owner_pdg_code);
        }
    }

    float global_true_drift_min = 1e5, global_true_drift_max = -1e5;
    float wire_min_u_truth = 1e5, wire_max_u_truth = -1e5;
    float wire_min_v_truth = 1e5, wire_max_v_truth = -1e5;
    float wire_min_w_truth = 1e5, wire_max_w_truth = -1e5;
    float buffer = 10.0;

    this->get_limits(true_hits_u_wire, true_hits_u_drift, wire_min_u_truth, wire_max_u_truth, global_true_drift_min, global_true_drift_max);
    this->get_limits(true_hits_v_wire, true_hits_v_drift, wire_min_v_truth, wire_max_v_truth, global_true_drift_min, global_true_drift_max);
    this->get_limits(true_hits_w_wire, true_hits_w_drift, wire_min_w_truth, wire_max_w_truth, global_true_drift_min, global_true_drift_max);

    TCanvas* canvas = new TCanvas("canvas", "", 1500, 1500);
    canvas->Divide(1, 3, 0, 0);

    TMultiGraph* mg_u = new TMultiGraph();
    TMultiGraph* mg_v = new TMultiGraph();
    TMultiGraph* mg_w = new TMultiGraph();

    mg_u->SetTitle(";Local Drift Coordinate;Local U Wire");
    mg_v->SetTitle(";Local Drift Coordinate;Local V Wire");
    mg_w->SetTitle(";Local Drift Coordinate;Local W Wire");

    /*TGraph* true_vertex_u = new TGraph();
    true_vertex_u->SetPoint(0, true_nu_vtx_x_, true_nu_vtx_u_wire_);
    true_vertex_u->SetMarkerStyle(4);
    true_vertex_u->SetMarkerColor(kBlack);
    true_vertex_u->SetMarkerSize(1.2);

    TGraph* true_vertex_v = new TGraph();
    true_vertex_v->SetPoint(0, true_nu_vtx_x_, true_nu_vtx_v_wire_);
    true_vertex_v->SetMarkerStyle(4);
    true_vertex_v->SetMarkerColor(kBlack);
    true_vertex_v->SetMarkerSize(1.2);

    TGraph* true_vertex_w = new TGraph();
    true_vertex_w->SetPoint(0, true_nu_vtx_x_, true_nu_vtx_w_wire_);
    true_vertex_w->SetMarkerStyle(4);
    true_vertex_w->SetMarkerColor(kBlack);
    true_vertex_w->SetMarkerSize(1.2);*/

    std::map<int, TGraph*> pdg_graphs;
    for (size_t i = 0; i < true_hits_u_wire.size(); ++i) {
        int pdg = std::abs(true_hits_u_owner.at(i));

        if (pdg_graphs.find(pdg) == pdg_graphs.end()) {
            pdg_graphs[pdg] = new TGraph();
            pdg_graphs[pdg]->SetMarkerStyle(20);
            pdg_graphs[pdg]->SetMarkerSize(0.5);
            pdg_graphs[pdg]->SetMarkerColor(kGray); // Default color

            if (pdg == 13) pdg_graphs[pdg]->SetMarkerColor(kBlue); // Muon
            else if (pdg == 11) pdg_graphs[pdg]->SetMarkerColor(kRed); // Electron
            else if (pdg == 2212) pdg_graphs[pdg]->SetMarkerColor(kGreen); // Proton
            else if (pdg == 211) pdg_graphs[pdg]->SetMarkerColor(kPink + 9); // Pion
            else if (pdg == 22) pdg_graphs[pdg]->SetMarkerColor(kOrange); // Photon
        }

        pdg_graphs[pdg]->SetPoint(pdg_graphs[pdg]->GetN(), true_hits_u_drift.at(i), true_hits_u_wire.at(i));
    }

    for (auto& entry : pdg_graphs) {
        mg_u->Add(entry.second);
    }
    //mg_u->Add(true_vertex_u);

    canvas->cd(1);
    mg_u->Draw("AP");
    mg_u->GetXaxis()->SetLimits(global_true_drift_min - buffer, global_true_drift_max + buffer);
    mg_u->GetXaxis()->SetTitleSize(0.05);  
    mg_u->GetYaxis()->SetTitleSize(0.05);

    pdg_graphs.clear(); 

    // Repeat for V view
    for (size_t i = 0; i < true_hits_v_wire.size(); ++i) {
        int pdg = std::abs(true_hits_v_owner.at(i));

        if (pdg_graphs.find(pdg) == pdg_graphs.end()) {
            pdg_graphs[pdg] = new TGraph();
            pdg_graphs[pdg]->SetMarkerStyle(20);
            pdg_graphs[pdg]->SetMarkerSize(0.5);
            pdg_graphs[pdg]->SetMarkerColor(kGray); // Default color

            if (pdg == 13) pdg_graphs[pdg]->SetMarkerColor(kBlue); // Muon
            else if (pdg == 11) pdg_graphs[pdg]->SetMarkerColor(kRed); // Electron
            else if (pdg == 2212) pdg_graphs[pdg]->SetMarkerColor(kGreen); // Proton
            else if (pdg == 211) pdg_graphs[pdg]->SetMarkerColor(kPink + 9); // Pion
            else if (pdg == 22) pdg_graphs[pdg]->SetMarkerColor(kOrange); // Photon
        }

        pdg_graphs[pdg]->SetPoint(pdg_graphs[pdg]->GetN(), true_hits_v_drift.at(i), true_hits_v_wire.at(i));
    }

    for (auto& entry : pdg_graphs) {
        mg_v->Add(entry.second);
    }
    //mg_v->Add(true_vertex_v);

    canvas->cd(2);
    mg_v->Draw("AP");
    mg_v->GetXaxis()->SetLimits(global_true_drift_min - buffer, global_true_drift_max + buffer);
    mg_v->GetXaxis()->SetTitleSize(0.05);  
    mg_v->GetYaxis()->SetTitleSize(0.05);

    pdg_graphs.clear(); 

    // Repeat for W view
    for (size_t i = 0; i < true_hits_w_wire.size(); ++i) {
        int pdg = std::abs(true_hits_w_owner.at(i));

        if (pdg_graphs.find(pdg) == pdg_graphs.end()) {
            pdg_graphs[pdg] = new TGraph();
            pdg_graphs[pdg]->SetMarkerStyle(20);
            pdg_graphs[pdg]->SetMarkerSize(0.5);
            pdg_graphs[pdg]->SetMarkerColor(kGray); // Default color

            if (pdg == 13) pdg_graphs[pdg]->SetMarkerColor(kBlue); // Muon
            else if (pdg == 11) pdg_graphs[pdg]->SetMarkerColor(kRed); // Electron
            else if (pdg == 2212) pdg_graphs[pdg]->SetMarkerColor(kGreen); // Proton
            else if (pdg == 211) pdg_graphs[pdg]->SetMarkerColor(kPink + 9); // Pion
            else if (pdg == 22) pdg_graphs[pdg]->SetMarkerColor(kOrange); // Photon
        }

        pdg_graphs[pdg]->SetPoint(pdg_graphs[pdg]->GetN(), true_hits_w_drift.at(i), true_hits_w_wire.at(i));
    }

    for (auto& entry : pdg_graphs) {
        mg_w->Add(entry.second);
    }
    //mg_w->Add(true_vertex_w);

    canvas->cd(3);
    mg_w->Draw("AP");
    mg_w->GetXaxis()->SetLimits(global_true_drift_min - buffer, global_true_drift_max + buffer);
    mg_w->GetXaxis()->SetTitleSize(0.05);  
    mg_w->GetYaxis()->SetTitleSize(0.05);
    canvas->SaveAs((filename + "_truth_hits.png").c_str());

    delete canvas;
    delete mg_u;
    delete mg_v;
    delete mg_w;
}

void VisualiseEventFilter::get_limits(const std::vector<float>& wire_coord_vec, const std::vector<float>& drift_coord_vec,
                   float& global_wire_min, float& global_wire_max, float& global_drift_min, float& global_drift_max) const
{
    if (!wire_coord_vec.empty() && !drift_coord_vec.empty()) {
        float local_wire_min = *std::min_element(wire_coord_vec.begin(), wire_coord_vec.end());
        float local_wire_max = *std::max_element(wire_coord_vec.begin(), wire_coord_vec.end());
        float local_drift_min = *std::min_element(drift_coord_vec.begin(), drift_coord_vec.end());
        float local_drift_max = *std::max_element(drift_coord_vec.begin(), drift_coord_vec.end());

        global_wire_min = std::min(global_wire_min, local_wire_min);
        global_wire_max = std::max(global_wire_max, local_wire_max);
        global_drift_min = std::min(global_drift_min, local_drift_min);
        global_drift_max = std::max(global_drift_max, local_drift_max);
    }
}

DEFINE_ART_MODULE(VisualiseEventFilter)
