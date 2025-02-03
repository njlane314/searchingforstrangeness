#ifndef VISUALISATION_H
#define VISUALISATION_H

#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "CommonFunctions/Pandora.h"

#include "SignatureTools/SignatureToolBase.h"
#include "SignatureTools/VertexToolBase.h"

#include "TVector3.h"
#include <vector>
#include <map>
#include <algorithm>
#include <string>

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

namespace common
{
    void visualiseTrueEvent(const art::Event& e,
                    const art::InputTag& mcp_producer,
                    const art::InputTag& hit_producer,
                    const art::InputTag& backtrack_tag,
                    const std::string& filename)
    {
        auto getLimits = [](const std::vector<float>& wire_coords, const std::vector<float>& drift_coords,
                    float& wire_min, float& wire_max, float& drift_min, float& drift_max)
        {
            if (!wire_coords.empty() && !drift_coords.empty()) {
                wire_min = std::min(wire_min, *std::min_element(wire_coords.begin(), wire_coords.end()));
                wire_max = std::max(wire_max, *std::max_element(wire_coords.begin(), wire_coords.end()));
                drift_min = std::min(drift_min, *std::min_element(drift_coords.begin(), drift_coords.end()));
                drift_max = std::max(drift_max, *std::max_element(drift_coords.begin(), drift_coords.end()));

                if ((wire_max - wire_min) < 100.0f) {
                    float padd = (100.0f - (wire_max - wire_min)) / 2.0f;
                    wire_min -= padd;
                    wire_max += padd;
                }

                if ((drift_max - drift_min) < 100.0f) {
                    float padd = (100.0f - (drift_max - drift_min)) / 2.0f;
                    drift_min -= padd;
                    drift_max += padd;
                }
            }
        };

        art::Handle<std::vector<simb::MCParticle>> mc_particle_handle; 
        std::vector<art::Ptr<simb::MCParticle>> mc_particle_vector;
        lar_pandora::MCParticleMap mc_particle_map;

        if (!e.getByLabel(mcp_producer, mc_particle_handle))
            throw cet::exception("Common") << "failed to find any mc particles in event" << std::endl;
        art::fill_ptr_vector(mc_particle_vector, mc_particle_handle);
        lar_pandora::LArPandoraHelper::BuildMCParticleMap(mc_particle_vector, mc_particle_map);

        art::Handle<std::vector<recob::Hit>> evt_hits;
        std::vector<art::Ptr<recob::Hit>> hit_vector;
        
        if (!e.getByLabel(hit_producer, evt_hits))
            throw cet::exception("Common") << "failed to find any hits in event" << std::endl;
        art::fill_ptr_vector(hit_vector, evt_hits);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(evt_hits, e, backtrack_tag);

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

        getLimits(true_hits_u_wire, true_hits_u_drift, wire_min_u_truth, wire_max_u_truth, global_true_drift_min, global_true_drift_max);
        getLimits(true_hits_v_wire, true_hits_v_drift, wire_min_v_truth, wire_max_v_truth, global_true_drift_min, global_true_drift_max);
        getLimits(true_hits_w_wire, true_hits_w_drift, wire_min_w_truth, wire_max_w_truth, global_true_drift_min, global_true_drift_max);

        TCanvas* canvas = new TCanvas("canvas", "", 1500, 1500);
        canvas->Divide(1, 3, 0, 0);

        TMultiGraph* mg_u = new TMultiGraph();
        TMultiGraph* mg_v = new TMultiGraph();
        TMultiGraph* mg_w = new TMultiGraph();

        mg_u->SetTitle(";Local Drift Coordinate;Local U Wire");
        mg_v->SetTitle(";Local Drift Coordinate;Local V Wire");
        mg_w->SetTitle(";Local Drift Coordinate;Local W Wire");

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
                else if (pdg == 321) pdg_graphs[pdg]->SetMarkerColor(kMagenta); // Kaon
                else if (pdg == 3222 || pdg == 3112) pdg_graphs[pdg]->SetMarkerColor(kCyan); // Sigma
            }

            pdg_graphs[pdg]->SetPoint(pdg_graphs[pdg]->GetN(), true_hits_u_drift.at(i), true_hits_u_wire.at(i));
        }

        for (auto& entry : pdg_graphs) {
            mg_u->Add(entry.second);
        }

        canvas->cd(1);
        mg_u->Draw("AP");
        mg_u->GetXaxis()->SetLimits(global_true_drift_min - buffer, global_true_drift_max + buffer);
        mg_u->GetYaxis()->SetRangeUser(wire_min_u_truth, wire_max_u_truth);
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
                else if (pdg == 321) pdg_graphs[pdg]->SetMarkerColor(kMagenta); // Kaon
                else if (pdg == 3222 || pdg == 3112) pdg_graphs[pdg]->SetMarkerColor(kCyan); // Sigma
            }

            pdg_graphs[pdg]->SetPoint(pdg_graphs[pdg]->GetN(), true_hits_v_drift.at(i), true_hits_v_wire.at(i));
        }

        for (auto& entry : pdg_graphs) {
            mg_v->Add(entry.second);
        }

        canvas->cd(2);
        mg_v->Draw("AP");
        mg_v->GetXaxis()->SetLimits(global_true_drift_min - buffer, global_true_drift_max + buffer);
        mg_v->GetYaxis()->SetRangeUser(wire_min_v_truth, wire_max_v_truth);
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
                else if (pdg == 321) pdg_graphs[pdg]->SetMarkerColor(kMagenta); // Kaon
                else if (pdg == 3222 || pdg == 3112) pdg_graphs[pdg]->SetMarkerColor(kCyan); // Sigma
            }

            pdg_graphs[pdg]->SetPoint(pdg_graphs[pdg]->GetN(), true_hits_w_drift.at(i), true_hits_w_wire.at(i));
        }

        for (auto& entry : pdg_graphs) {
            mg_w->Add(entry.second);
        }

        canvas->cd(3);
        mg_w->Draw("AP");
        mg_w->GetXaxis()->SetLimits(global_true_drift_min - buffer, global_true_drift_max + buffer);
        mg_w->GetYaxis()->SetRangeUser(wire_min_w_truth, wire_max_w_truth);
        mg_w->GetXaxis()->SetTitleSize(0.05);  
        mg_w->GetYaxis()->SetTitleSize(0.05);
        canvas->SaveAs((filename + "_truth_hits.png").c_str());

        delete canvas;
        delete mg_u;
        delete mg_v;
        delete mg_w;
    }

    /*void visualisePandoraEvent()
    {
        auto getLimits = [](const std::vector<float>& wire_coords, const std::vector<float>& drift_coords,
                    float& wire_min, float& wire_max, float& drift_min, float& drift_max)
        {
            if (!wire_coords.empty() && !drift_coords.empty()) {
                wire_min = std::min(wire_min, *std::min_element(wire_coords.begin(), wire_coords.end()));
                wire_max = std::max(wire_max, *std::max_element(wire_coords.begin(), wire_coords.end()));
                drift_min = std::min(drift_min, *std::min_element(drift_coords.begin(), drift_coords.end()));
                drift_max = std::max(drift_max, *std::max_element(drift_coords.begin(), drift_coords.end()));

                if ((wire_max - wire_min) < 100.0f) {
                    float padd = (100.0f - (wire_max - wire_min)) / 2.0f;
                    wire_min -= padd;
                    wire_max += padd;
                }

                if ((drift_max - drift_min) < 100.0f) {
                    float padd = (100.0f - (drift_max - drift_min)) / 2.0f;
                    drift_min -= padd;
                    drift_max += padd;
                }
            }
        };

        float global_reco_drift_min = 1e10, global_reco_drift_max = -1e10;
        float wire_min_u_slice = 1e10, wire_max_u_slice = -1e10;
        float wire_min_v_slice = 1e10, wire_max_v_slice = -1e10;
        float wire_min_w_slice = 1e10, wire_max_w_slice = -1e10;
        float buffer = 10.0;

        // Calculate global min and max for reconstructed hits
        for (size_t i = 0; i < reco_hits_u_wire_->size(); ++i) {
            get_limits(reco_hits_u_wire_->at(i), reco_hits_u_drift_->at(i), wire_min_u_slice, wire_max_u_slice, global_reco_drift_min, global_reco_drift_max);
        }
        for (size_t i = 0; i < reco_hits_v_wire_->size(); ++i) {
            get_limits(reco_hits_v_wire_->at(i), reco_hits_v_drift_->at(i), wire_min_v_slice, wire_max_v_slice, global_reco_drift_min, global_reco_drift_max);
        }
        for (size_t i = 0; i < reco_hits_w_wire_->size(); ++i) {
            get_limits(reco_hits_w_wire_->at(i), reco_hits_w_drift_->at(i), wire_min_w_slice, wire_max_w_slice, global_reco_drift_min, global_reco_drift_max);
        }

        // Create TMultiGraphs and TGraphs for each view (U, V, W)
        TCanvas* c5 = new TCanvas("c5", "", 1500, 1500);
        c5->Divide(1, 3, 0, 0);

        TMultiGraph* reco_mg_u = new TMultiGraph();
        TMultiGraph* reco_mg_v = new TMultiGraph();
        TMultiGraph* reco_mg_w = new TMultiGraph();

        reco_mg_u->SetTitle(";Local Drift Coordinate;Local U Wire");
        reco_mg_v->SetTitle(";Local Drift Coordinate;Local V Wire");
        reco_mg_w->SetTitle(";Local Drift Coordinate;Local W Wire");

        std::vector<int> color_map = {
            kMagenta, kCyan, kYellow, kAzure, kSpring, kTeal, kRose, kGray, kBlack, kViolet,
            kOrange + 7, kBlue - 9, kGreen + 3, kViolet + 9, kCyan + 3, kYellow + 2, kGray + 2
        };

        int colour_map_index = 0;
        for (size_t i = 0; i < reco_hits_u_drift_->size(); ++i) {
            int particle_color = color_map[colour_map_index];
            colour_map_index++;

            TGraph* pfp_graph_u = new TGraph();
            pfp_graph_u->SetMarkerStyle(20);
            pfp_graph_u->SetMarkerSize(0.5);
            pfp_graph_u->SetMarkerColor(particle_color);

            for (size_t hit = 0; hit < reco_hits_u_drift_->at(i).size(); ++hit) {
                pfp_graph_u->SetPoint(
                    pfp_graph_u->GetN(), 
                    reco_hits_u_drift_->at(i).at(hit), 
                    reco_hits_u_wire_->at(i).at(hit)
                );
            }

            reco_mg_u->Add(pfp_graph_u);
        }

        colour_map_index = 0;
        for (size_t i = 0; i < reco_hits_v_drift_->size(); ++i) {
            int particle_color = color_map[colour_map_index];
            colour_map_index++;

            TGraph* pfp_graph_v = new TGraph();
            pfp_graph_v->SetMarkerStyle(20);
            pfp_graph_v->SetMarkerSize(0.5);
            pfp_graph_v->SetMarkerColor(particle_color);

            for (size_t hit = 0; hit < reco_hits_v_drift_->at(i).size(); ++hit) {
                pfp_graph_v->SetPoint(
                    pfp_graph_v->GetN(), 
                    reco_hits_v_drift_->at(i).at(hit), 
                    reco_hits_v_wire_->at(i).at(hit)
                );
            }

            reco_mg_v->Add(pfp_graph_v);
        }

        colour_map_index = 0;
        for (size_t i = 0; i < reco_hits_w_drift_->size(); ++i) {
            int particle_color = color_map[colour_map_index];
            colour_map_index++;

            TGraph* pfp_graph_w = new TGraph();
            pfp_graph_w->SetMarkerStyle(20);
            pfp_graph_w->SetMarkerSize(0.5);
            pfp_graph_w->SetMarkerColor(particle_color);

            for (size_t hit = 0; hit < reco_hits_w_drift_->at(i).size(); ++hit) {
                pfp_graph_w->SetPoint(pfp_graph_w->GetN(), reco_hits_w_drift_->at(i).at(hit), reco_hits_w_wire_->at(i).at(hit));
            }

            reco_mg_w->Add(pfp_graph_w);
        }

        c5->cd(1);
        reco_mg_u->Draw("AP");
        reco_mg_u->GetXaxis()->SetLimits(global_reco_drift_min - buffer, global_reco_drift_max + buffer);
        reco_mg_u->GetXaxis()->SetTitleSize(0.05);  
        reco_mg_u->GetYaxis()->SetTitleSize(0.05);

        c5->cd(2);
        reco_mg_v->Draw("AP");
        reco_mg_v->GetXaxis()->SetLimits(global_reco_drift_min - buffer, global_reco_drift_max + buffer);
        reco_mg_v->GetXaxis()->SetTitleSize(0.05);  
        reco_mg_v->GetYaxis()->SetTitleSize(0.05);

        c5->cd(3);
        reco_mg_w->Draw("AP");
        reco_mg_w->GetXaxis()->SetLimits(global_reco_drift_min - buffer, global_reco_drift_max + buffer);
        reco_mg_w->GetXaxis()->SetTitleSize(0.05);  
        reco_mg_w->GetYaxis()->SetTitleSize(0.05);

        std::string filename = "reco_interaction_hits_" + std::to_string(run_) + "_" + std::to_string(subrun_) + "_" + std::to_string(event_);
        c5->SaveAs(("./plots/" + filename + ".pdf").c_str());
    }*/

    /*void visualiseSignature(const art::Event& e,
                        const art::InputTag& mcp_producer,
                        const art::InputTag& hit_producer,
                        const art::InputTag& backtrack_tag,
                        const signature::Pattern& patt,
                        const std::string& filename)
    {
        auto getLimits = [](const std::vector<float>& wire_coords, const std::vector<float>& drift_coords,
                    float& wire_min, float& wire_max, float& drift_min, float& drift_max)
        {
            if (!wire_coords.empty() && !drift_coords.empty()) {
                wire_min = std::min(wire_min, *std::min_element(wire_coords.begin(), wire_coords.end()));
                wire_max = std::max(wire_max, *std::max_element(wire_coords.begin(), wire_coords.end()));
                drift_min = std::min(drift_min, *std::min_element(drift_coords.begin(), drift_coords.end()));
                drift_max = std::max(drift_max, *std::max_element(drift_coords.begin(), drift_coords.end()));

                if ((wire_max - wire_min) < 100.0f) {
                    float padd = (100.0f - (wire_max - wire_min)) / 2.0f;
                    wire_min -= padd;
                    wire_max += padd;
                }

                if ((drift_max - drift_min) < 100.0f) {
                    float padd = (100.0f - (drift_max - drift_min)) / 2.0f;
                    drift_min -= padd;
                    drift_max += padd;
                }
            }
        };

        art::Handle<std::vector<simb::MCParticle>> mc_particle_handle; 
        std::vector<art::Ptr<simb::MCParticle>> mc_particle_vector;
        lar_pandora::MCParticleMap mc_particle_map;

        if (!e.getByLabel(mcp_producer, mc_particle_handle))
            throw cet::exception("Common") << "failed to find any mc particles in event" << std::endl;
        art::fill_ptr_vector(mc_particle_vector, mc_particle_handle);
        lar_pandora::LArPandoraHelper::BuildMCParticleMap(mc_particle_vector, mc_particle_map);

        art::Handle<std::vector<recob::Hit>> evt_hits;
        std::vector<art::Ptr<recob::Hit>> hit_vector;
        
        if (!e.getByLabel(hit_producer, evt_hits))
            throw cet::exception("Common") << "failed to find any hits in event" << std::endl;
        art::fill_ptr_vector(hit_vector, evt_hits);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(evt_hits, e, backtrack_tag);

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

        std::vector<float> sig_u_drift, sig_u_wire, other_u_drift, other_u_wire;
        std::vector<float> sig_v_drift, sig_v_wire, other_v_drift, other_v_wire;
        std::vector<float> sig_w_drift, sig_w_wire, other_w_drift, other_w_wire;

        std::unordered_set<int> signature_track_ids;
        for (const auto& signature : patt) {
            for (const auto& mcp : signature) {
                signature_track_ids.insert(mcp->TrackId());
            }
        }

        for (const art::Ptr<recob::Hit> &hit : hit_vector)
        {
            bool is_signature_hit = false;

            const auto& mcps = assoc_mc_part.at(hit.key());
            const auto& match_data = assoc_mc_part.data(hit.key());

            auto hit_to_track_it = hits_to_track_map.find(hit.key());
            if (hit_to_track_it == hits_to_track_map.end()) {
                continue; 
            }

            for (size_t i = 0; i < mcps.size(); ++i) {
                if (signature_track_ids.count(mcps[i]->TrackId())) {
                    is_signature_hit = true;
                    break;
                }
            }

            common::PandoraView pandora_view = common::GetPandoraView(hit);
            TVector3 pandora_pos = common::GetPandoraHitPosition(e, hit, pandora_view);

            if (pandora_view == common::TPC_VIEW_U) {
                (is_signature_hit ? sig_u_drift : other_u_drift).push_back(pandora_pos.X());
                (is_signature_hit ? sig_u_wire : other_u_wire).push_back(pandora_pos.Z());
            }
            else if (pandora_view == common::TPC_VIEW_V) {
                (is_signature_hit ? sig_v_drift : other_v_drift).push_back(pandora_pos.X());
                (is_signature_hit ? sig_v_wire : other_v_wire).push_back(pandora_pos.Z());
            }
            else if (pandora_view == common::TPC_VIEW_W) {
                (is_signature_hit ? sig_w_drift : other_w_drift).push_back(pandora_pos.X());
                (is_signature_hit ? sig_w_wire : other_w_wire).push_back(pandora_pos.Z());
            }
        }

        float drift_min = 1e5, drift_max = -1e5;
        float wire_min_u = 1e5, wire_max_u = -1e5;
        float wire_min_v = 1e5, wire_max_v = -1e5;
        float wire_min_w = 1e5, wire_max_w = -1e5;
        float buffer = 10.0;

        getLimits(sig_u_wire, sig_u_drift, wire_min_u, wire_max_u, drift_min, drift_max);
        getLimits(other_u_wire, other_u_drift, wire_min_u, wire_max_u, drift_min, drift_max);
        getLimits(sig_v_wire, sig_v_drift, wire_min_v, wire_max_v, drift_min, drift_max);
        getLimits(other_v_wire, other_v_drift, wire_min_v, wire_max_v, drift_min, drift_max);
        getLimits(sig_w_wire, sig_w_drift, wire_min_w, wire_max_w, drift_min, drift_max);
        getLimits(other_w_wire, other_w_drift, wire_min_w, wire_max_w, drift_min, drift_max);

        TCanvas* canvas = new TCanvas("canvas", "", 1500, 1500);
        canvas->Divide(1, 3, 0, 0);

        TMultiGraph* mg_u = new TMultiGraph();
        TMultiGraph* mg_v = new TMultiGraph();
        TMultiGraph* mg_w = new TMultiGraph();

        mg_u->SetTitle(";Local Drift Coordinate;Local U Wire");
        mg_v->SetTitle(";Local Drift Coordinate;Local V Wire");
        mg_w->SetTitle(";Local Drift Coordinate;Local W Wire");

        TGraph* sig_u = new TGraph();
        sig_u->SetMarkerStyle(20);
        sig_u->SetMarkerSize(0.5);
        sig_u->SetMarkerColor(kGreen);
        for (size_t i = 0; i < sig_u_wire.size(); ++i) 
            sig_u->SetPoint(sig_u->GetN(), sig_u_drift.at(i), sig_u_wire.at(i));

        TGraph* other_u = new TGraph();
        other_u->SetMarkerStyle(20);
        other_u->SetMarkerSize(0.5);
        other_u->SetMarkerColor(kGray);
        for (size_t i = 0; i < other_u_wire.size(); ++i) 
            other_u->SetPoint(other_u->GetN(), other_u_drift.at(i), other_u_wire.at(i));

        mg_u->Add(sig_u);
        mg_u->Add(other_u);

        TGraph* sig_v = new TGraph();
        sig_v->SetMarkerStyle(20);
        sig_v->SetMarkerSize(0.5);
        sig_v->SetMarkerColor(kGreen);
        for (size_t i = 0; i < sig_v_wire.size(); ++i) 
            sig_v->SetPoint(sig_v->GetN(), sig_v_drift.at(i), sig_v_wire.at(i));

        TGraph* other_v = new TGraph();
        other_v->SetMarkerStyle(20);
        other_v->SetMarkerSize(0.5); 
        other_v->SetMarkerColor(kGray);
        for (size_t i = 0; i < other_v_wire.size(); ++i) 
            other_v->SetPoint(other_v->GetN(), other_v_drift.at(i), other_v_wire.at(i));

        mg_v->Add(sig_v);
        mg_v->Add(other_v);

        TGraph* sig_w = new TGraph();
        sig_w->SetMarkerStyle(20);
        sig_w->SetMarkerSize(0.5);
        sig_w->SetMarkerColor(kGreen);
        for (size_t i = 0; i < sig_w_wire.size(); ++i) 
            sig_w->SetPoint(sig_w->GetN(), sig_w_drift.at(i), sig_w_wire.at(i));

        TGraph* other_w = new TGraph();
        other_w->SetMarkerStyle(20);
        other_w->SetMarkerSize(0.5);
        other_w->SetMarkerColor(kGray);
        for (size_t i = 0; i < other_w_wire.size(); ++i) 
            other_w->SetPoint(other_w->GetN(), other_w_drift.at(i), other_w_wire.at(i));

        mg_w->Add(sig_w);
        mg_w->Add(other_w);

        canvas->cd(1);
        mg_u->Draw("AP");
        mg_u->GetXaxis()->SetLimits(drift_min - buffer, drift_max + buffer);
        mg_u->GetYaxis()->SetRangeUser(wire_min_u, wire_max_u);
        mg_u->GetXaxis()->SetTitleSize(0.05);  
        mg_u->GetYaxis()->SetTitleSize(0.05);

        canvas->cd(2);
        mg_v->Draw("AP");
        mg_v->GetXaxis()->SetLimits(drift_min - buffer, drift_max + buffer);
        mg_v->GetYaxis()->SetRangeUser(wire_min_v, wire_max_v);
        mg_v->GetXaxis()->SetTitleSize(0.05);  
        mg_v->GetYaxis()->SetTitleSize(0.05);

        canvas->cd(3);
        mg_w->Draw("AP");
        mg_w->GetXaxis()->SetLimits(drift_min - buffer, drift_max + buffer);
        mg_w->GetYaxis()->SetRangeUser(wire_min_w, wire_max_w);
        mg_w->GetXaxis()->SetTitleSize(0.05);  
        mg_w->GetYaxis()->SetTitleSize(0.05);

        canvas->SaveAs((filename + "_signature_hits.png").c_str());

        delete canvas;
        delete mg_u;
        delete mg_v;
        delete mg_w;
    }*/


    void visualiseSignature(const art::Event& e,
                    const art::InputTag& mcp_producer,
                    const art::InputTag& hit_producer,
                    const art::InputTag& backtrack_tag,
                    const signature::Pattern& patt,
                    const std::string& filename)
    {
        auto getLimits = [](const std::vector<float>& wire_coords, const std::vector<float>& drift_coords,
                    float& wire_min, float& wire_max, float& drift_min, float& drift_max)
        {
            if (!wire_coords.empty() && !drift_coords.empty()) {
                wire_min = std::min(wire_min, *std::min_element(wire_coords.begin(), wire_coords.end()));
                wire_max = std::max(wire_max, *std::max_element(wire_coords.begin(), wire_coords.end()));
                drift_min = std::min(drift_min, *std::min_element(drift_coords.begin(), drift_coords.end()));
                drift_max = std::max(drift_max, *std::max_element(drift_coords.begin(), drift_coords.end()));

                if ((wire_max - wire_min) < 100.0f) {
                    float padd = (100.0f - (wire_max - wire_min)) / 2.0f;
                    wire_min -= padd;
                    wire_max += padd;
                }

                if ((drift_max - drift_min) < 100.0f) {
                    float padd = (100.0f - (drift_max - drift_min)) / 2.0f;
                    drift_min -= padd;
                    drift_max += padd;
                }
            }
        };

        art::Handle<std::vector<simb::MCParticle>> mc_particle_handle; 
        std::vector<art::Ptr<simb::MCParticle>> mc_particle_vector;
        lar_pandora::MCParticleMap mc_particle_map;

        if (!e.getByLabel(mcp_producer, mc_particle_handle))
            throw cet::exception("Common") << "failed to find any mc particles in event" << std::endl;
        art::fill_ptr_vector(mc_particle_vector, mc_particle_handle);
        lar_pandora::LArPandoraHelper::BuildMCParticleMap(mc_particle_vector, mc_particle_map);

        art::Handle<std::vector<recob::Hit>> evt_hits;
        std::vector<art::Ptr<recob::Hit>> hit_vector;
        
        if (!e.getByLabel(hit_producer, evt_hits))
            throw cet::exception("Common") << "failed to find any hits in event" << std::endl;
        art::fill_ptr_vector(hit_vector, evt_hits);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assoc_mc_part = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(evt_hits, e, backtrack_tag);

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

            int owner_trackid = mc_particle_map.at(hit_to_track_it->second)->TrackId();

            if (pandora_view == common::TPC_VIEW_U) {
                true_hits_u_wire.push_back(pandora_pos.Z());
                true_hits_u_drift.push_back(pandora_pos.X());
                true_hits_u_owner.push_back(owner_trackid);
            }
            else if (pandora_view == common::TPC_VIEW_V) {
                true_hits_v_wire.push_back(pandora_pos.Z());
                true_hits_v_drift.push_back(pandora_pos.X());
                true_hits_v_owner.push_back(owner_trackid);
            }
            else if (pandora_view == common::TPC_VIEW_W) {
                true_hits_w_wire.push_back(pandora_pos.Z());
                true_hits_w_drift.push_back(pandora_pos.X());
                true_hits_w_owner.push_back(owner_trackid);
            }
        }

        float global_true_drift_min = 1e5, global_true_drift_max = -1e5;
        float wire_min_u_truth = 1e5, wire_max_u_truth = -1e5;
        float wire_min_v_truth = 1e5, wire_max_v_truth = -1e5;
        float wire_min_w_truth = 1e5, wire_max_w_truth = -1e5;
        float buffer = 10.0;

        getLimits(true_hits_u_wire, true_hits_u_drift, wire_min_u_truth, wire_max_u_truth, global_true_drift_min, global_true_drift_max);
        getLimits(true_hits_v_wire, true_hits_v_drift, wire_min_v_truth, wire_max_v_truth, global_true_drift_min, global_true_drift_max);
        getLimits(true_hits_w_wire, true_hits_w_drift, wire_min_w_truth, wire_max_w_truth, global_true_drift_min, global_true_drift_max);

        TCanvas* canvas = new TCanvas("canvas", "", 1500, 1500);
        canvas->Divide(1, 3, 0, 0);

        TMultiGraph* mg_u = new TMultiGraph();
        TMultiGraph* mg_v = new TMultiGraph();
        TMultiGraph* mg_w = new TMultiGraph();

        mg_u->SetTitle(";Local Drift Coordinate;Local U Wire");
        mg_v->SetTitle(";Local Drift Coordinate;Local V Wire");
        mg_w->SetTitle(";Local Drift Coordinate;Local W Wire");

        TGraph* sig_u = new TGraph();
        sig_u->SetMarkerStyle(20);
        sig_u->SetMarkerSize(0.5);
        sig_u->SetMarkerColor(kGreen);

        TGraph* back_u = new TGraph();
        back_u->SetMarkerStyle(20);
        back_u->SetMarkerSize(0.5);
        back_u->SetMarkerColor(kGray);
        for (size_t i = 0; i < true_hits_u_wire.size(); ++i) {
            int trackid = std::abs(true_hits_u_owner.at(i));

            bool is_sig = false;
            for (const auto& signature : patt) {
                for (const auto& mcp : signature.second) {
                    if (mcp->TrackId() == trackid){
                        is_sig = true; 
                    }
                }
            }

            if (is_sig)
                sig_u->SetPoint(sig_u->GetN(), true_hits_u_drift.at(i), true_hits_u_wire.at(i));
            else if (!is_sig) 
                back_u->SetPoint(back_u->GetN(), true_hits_u_drift.at(i), true_hits_u_wire.at(i));
        }

        mg_u->Add(sig_u); 
        mg_u->Add(back_u);

        canvas->cd(1);
        mg_u->Draw("AP");
        mg_u->GetXaxis()->SetLimits(global_true_drift_min - buffer, global_true_drift_max + buffer);
        mg_u->GetYaxis()->SetRangeUser(wire_min_u_truth, wire_max_u_truth);
        mg_u->GetXaxis()->SetTitleSize(0.05);  
        mg_u->GetYaxis()->SetTitleSize(0.05);

        TGraph* sig_v = new TGraph();
        sig_v->SetMarkerStyle(20);
        sig_v->SetMarkerSize(0.5);
        sig_v->SetMarkerColor(kGreen);

        TGraph* back_v = new TGraph();
        back_v->SetMarkerStyle(20);
        back_v->SetMarkerSize(0.5);
        back_v->SetMarkerColor(kGray);
        for (size_t i = 0; i < true_hits_v_wire.size(); ++i) {
            int trackid = std::abs(true_hits_v_owner.at(i));

            bool is_sig = false;
            for (const auto& signature : patt) {
                for (const auto& mcp : signature.second) {
                    if (mcp->TrackId() == trackid){
                        is_sig = true; 
                    }
                }
            }

            if (is_sig)
                sig_v->SetPoint(sig_v->GetN(), true_hits_v_drift.at(i), true_hits_v_wire.at(i));
            else if (!is_sig) 
                back_v->SetPoint(back_v->GetN(), true_hits_v_drift.at(i), true_hits_v_wire.at(i));
        }

        mg_v->Add(sig_v); 
        mg_v->Add(back_v);

        canvas->cd(2);
        mg_v->Draw("AP");
        mg_v->GetXaxis()->SetLimits(global_true_drift_min - buffer, global_true_drift_max + buffer);
        mg_v->GetYaxis()->SetRangeUser(wire_min_v_truth, wire_max_v_truth);
        mg_v->GetXaxis()->SetTitleSize(0.05);  
        mg_v->GetYaxis()->SetTitleSize(0.05);


        TGraph* sig_w = new TGraph();
        sig_w->SetMarkerStyle(20);
        sig_w->SetMarkerSize(0.5);
        sig_w->SetMarkerColor(kGreen);

        TGraph* back_w = new TGraph();
        back_w->SetMarkerStyle(20);
        back_w->SetMarkerSize(0.5);
        back_w->SetMarkerColor(kGray);
        for (size_t i = 0; i < true_hits_w_wire.size(); ++i) {
            int trackid = std::abs(true_hits_w_owner.at(i));

            bool is_sig = false;
            for (const auto& signature : patt) {
                for (const auto& mcp : signature.second) {
                    if (mcp->TrackId() == trackid){
                        is_sig = true; 
                    }
                }
            }

            if (is_sig)
                sig_w->SetPoint(sig_w->GetN(), true_hits_w_drift.at(i), true_hits_w_wire.at(i));
            else if (!is_sig) 
                back_w->SetPoint(back_w->GetN(), true_hits_w_drift.at(i), true_hits_w_wire.at(i));
        }

        mg_w->Add(sig_w); 
        mg_w->Add(back_w);

        canvas->cd(3);
        mg_w->Draw("AP");
        mg_w->GetXaxis()->SetLimits(global_true_drift_min - buffer, global_true_drift_max + buffer);
        mg_w->GetYaxis()->SetRangeUser(wire_min_w_truth, wire_max_w_truth);
        mg_w->GetXaxis()->SetTitleSize(0.05);  
        mg_w->GetYaxis()->SetTitleSize(0.05);
        canvas->SaveAs((filename + "_signature_hits.png").c_str());

        delete canvas;
        delete mg_u;
        delete mg_v;
        delete mg_w;
    }

}

#endif
