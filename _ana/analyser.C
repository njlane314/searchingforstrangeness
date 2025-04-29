#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

void analyse(const std::vector<std::string>& in_file_name, const std::string& output_filename) {
    TChain fiducial_ch("fiducialFilter/EventSelectionFilter");
    TChain topology_ch("topologyTrainingFilter/EventSelectionFilter");
    for (const auto& f_name : in_file_name) {
        fiducial_ch.Add(f_name.c_str());
        topology_ch.Add(f_name.c_str());
    }

    const Long64_t multiplier = 1000000; // Adjust based on your data
    std::string minor_expr = "sub * " + std::to_string(multiplier) + " + evt";
    fiducial_ch.BuildIndex("run", minor_expr.c_str());
    topology_ch.BuildIndex("run", minor_expr.c_str());

    int temp_run, temp_sub, temp_evt;
    std::set<std::tuple<int, int, int>> unique_events;
    for (auto& chain : {&fiducial_ch, &topology_ch}) {
        chain->SetBranchAddress("run", &temp_run);
        chain->SetBranchAddress("sub", &temp_sub);
        chain->SetBranchAddress("evt", &temp_evt);
        for (Long64_t i = 0; i < chain->GetEntries(); ++i) {
            chain->GetEntry(i);
            unique_events.insert(std::make_tuple(temp_run, temp_sub, temp_evt));
        }
    }

    int fid_run, fid_sub, fid_evt, fid_selected;
    double fid_reco_nu_vtx_x, fid_reco_nu_vtx_y, fid_reco_nu_vtx_z;
    int fid_event_category;
    float fid_mcf_nu_e, fid_mcf_lep_e;
    int fid_mcf_actvol, fid_mcf_fidvol, fid_mcf_nmm, fid_mcf_nmp, fid_mcf_nem, fid_mcf_nep, fid_mcf_np0, fid_mcf_npp, fid_mcf_npm, fid_mcf_nkp, fid_mcf_nkm, fid_mcf_nk0, fid_mcf_npr, fid_mcf_nne, fid_mcf_nlambda, fid_mcf_nsigma_p, fid_mcf_nsigma_0, fid_mcf_nsigma_m, fid_mcf_nxi_0, fid_mcf_nxi_m, fid_mcf_nomega;

    fiducial_ch.SetBranchAddress("run", &fid_run);
    fiducial_ch.SetBranchAddress("sub", &fid_sub);
    fiducial_ch.SetBranchAddress("evt", &fid_evt);
    fiducial_ch.SetBranchAddress("selected", &fid_selected);
    fiducial_ch.SetBranchAddress("reco_nu_vtx_x", &fid_reco_nu_vtx_x);
    fiducial_ch.SetBranchAddress("reco_nu_vtx_y", &fid_reco_nu_vtx_y);
    fiducial_ch.SetBranchAddress("reco_nu_vtx_z", &fid_reco_nu_vtx_z);
    fiducial_ch.SetBranchAddress("event_category", &fid_event_category);
    fiducial_ch.SetBranchAddress("mcf_nu_e", &fid_mcf_nu_e);
    fiducial_ch.SetBranchAddress("mcf_lep_e", &fid_mcf_lep_e);
    fiducial_ch.SetBranchAddress("mcf_actvol", &fid_mcf_actvol);
    fiducial_ch.SetBranchAddress("mcf_fidvol", &fid_mcf_fidvol);
    fiducial_ch.SetBranchAddress("mcf_nmm", &fid_mcf_nmm);
    fiducial_ch.SetBranchAddress("mcf_nmp", &fid_mcf_nmp);
    fiducial_ch.SetBranchAddress("mcf_nem", &fid_mcf_nem);
    fiducial_ch.SetBranchAddress("mcf_nep", &fid_mcf_nep);
    fiducial_ch.SetBranchAddress("mcf_np0", &fid_mcf_np0);
    fiducial_ch.SetBranchAddress("mcf_npp", &fid_mcf_npp);
    fiducial_ch.SetBranchAddress("mcf_npm", &fid_mcf_npm);
    fiducial_ch.SetBranchAddress("mcf_nkp", &fid_mcf_nkp);
    fiducial_ch.SetBranchAddress("mcf_nkm", &fid_mcf_nkm);
    fiducial_ch.SetBranchAddress("mcf_nk0", &fid_mcf_nk0);
    fiducial_ch.SetBranchAddress("mcf_npr", &fid_mcf_npr);
    fiducial_ch.SetBranchAddress("mcf_nne", &fid_mcf_nne);
    fiducial_ch.SetBranchAddress("mcf_nlambda", &fid_mcf_nlambda);
    fiducial_ch.SetBranchAddress("mcf_nsigma_p", &fid_mcf_nsigma_p);
    fiducial_ch.SetBranchAddress("mcf_nsigma_0", &fid_mcf_nsigma_0);
    fiducial_ch.SetBranchAddress("mcf_nsigma_m", &fid_mcf_nsigma_m);
    fiducial_ch.SetBranchAddress("mcf_nxi_0", &fid_mcf_nxi_0);
    fiducial_ch.SetBranchAddress("mcf_nxi_m", &fid_mcf_nxi_m);
    fiducial_ch.SetBranchAddress("mcf_nomega", &fid_mcf_nomega);

    int top_run, top_sub, top_evt, top_selected;
    int top_nhits_u, top_nhit_v, top_nhit_w;
    float top_charge_u, top_charge_v, top_charge_w;
    float top_wirerange_u, top_wirerange_v, top_wirerange_w;
    float top_timerange_u, top_timerange_v, top_timerange_w;
    int top_event_category;
    float top_mcf_nu_e, top_mcf_lep_e;
    int top_mcf_actvol, top_mcf_fidvol, top_mcf_nmm, top_mcf_nmp, top_mcf_nem, top_mcf_nep, top_mcf_np0, top_mcf_npp, top_mcf_npm, top_mcf_nkp, top_mcf_nkm, top_mcf_nk0, top_mcf_npr, top_mcf_nne, top_mcf_nlambda, top_mcf_nsigma_p, top_mcf_nsigma_0, top_mcf_nsigma_m, top_mcf_nxi_0, top_mcf_nxi_m, top_mcf_nomega;

    topology_ch.SetBranchAddress("run", &top_run);
    topology_ch.SetBranchAddress("sub", &top_sub);
    topology_ch.SetBranchAddress("evt", &top_evt);
    topology_ch.SetBranchAddress("selected", &top_selected);
    topology_ch.SetBranchAddress("nhits_u", &top_nhits_u);
    topology_ch.SetBranchAddress("nhit_v", &top_nhit_v);
    topology_ch.SetBranchAddress("nhit_w", &top_nhit_w);
    topology_ch.SetBranchAddress("charge_u", &top_charge_u);
    topology_ch.SetBranchAddress("charge_v", &top_charge_v);
    topology_ch.SetBranchAddress("chare_w", &top_charge_w);
    topology_ch.SetBranchAddress("wirerange_u", &top_wirerange_u);
    topology_ch.SetBranchAddress("wirerange_v", &top_wirerange_v);
    topology_ch.SetBranchAddress("wirerange_w", &top_wirerange_w);
    topology_ch.SetBranchAddress("timerange_u", &top_timerange_u);
    topology_ch.SetBranchAddress("timerange_v", &top_timerange_v);
    topology_ch.SetBranchAddress("timerange_w", &top_timerange_w);
    topology_ch.SetBranchAddress("event_category", &top_event_category);
    topology_ch.SetBranchAddress("mcf_nu_e", &top_mcf_nu_e);
    topology_ch.SetBranchAddress("mcf_lep_e", &top_mcf_lep_e);
    topology_ch.SetBranchAddress("mcf_actvol", &top_mcf_actvol);
    topology_ch.SetBranchAddress("mcf_fidvol", &top_mcf_fidvol);
    topology_ch.SetBranchAddress("mcf_nmm", &top_mcf_nmm);
    topology_ch.SetBranchAddress("mcf_nmp", &top_mcf_nmp);
    topology_ch.SetBranchAddress("mcf_nem", &top_mcf_nem);
    topology_ch.SetBranchAddress("mcf_nep", &top_mcf_nep);
    topology_ch.SetBranchAddress("mcf_np0", &top_mcf_np0);
    topology_ch.SetBranchAddress("mcf_npp", &top_mcf_npp);
    topology_ch.SetBranchAddress("mcf_npm", &top_mcf_npm);
    topology_ch.SetBranchAddress("mcf_nkp", &top_mcf_nkp);
    topology_ch.SetBranchAddress("mcf_nkm", &top_mcf_nkm);
    topology_ch.SetBranchAddress("mcf_nk0", &top_mcf_nk0);
    topology_ch.SetBranchAddress("mcf_npr", &top_mcf_npr);
    topology_ch.SetBranchAddress("mcf_nne", &top_mcf_nne);
    topology_ch.SetBranchAddress("mcf_nlambda", &top_mcf_nlambda);
    topology_ch.SetBranchAddress("mcf_nsigma_p", &top_mcf_nsigma_p);
    topology_ch.SetBranchAddress("mcf_nsigma_0", &top_mcf_nsigma_0);
    topology_ch.SetBranchAddress("mcf_nsigma_m", &top_mcf_nsigma_m);
    topology_ch.SetBranchAddress("mcf_nxi_0", &top_mcf_nxi_0);
    topology_ch.SetBranchAddress("mcf_nxi_m", &top_mcf_nxi_m);
    topology_ch.SetBranchAddress("mcf_nomega", &top_mcf_nomega);

    TFile* out_file = new TFile(output_filename.c_str(), "recreate");
    TTree* out_tree = new TTree("strange_tree", "Strange Analysis Tree");

    int out_run, out_sub, out_evt;
    int out_selected_fiducial, out_selected_topology;
    double out_reco_nu_vtx_x, out_reco_nu_vtx_y, out_reco_nu_vtx_z;
    int out_nhits_u, out_nhit_v, out_nhit_w;
    float out_charge_u, out_charge_v, out_charge_w;
    float out_wirerange_u, out_wirerange_v, out_wirerange_w;
    float out_timerange_u, out_timerange_v, out_timerange_w;
    int out_event_category;
    float out_mcf_nu_e, out_mcf_lep_e;
    int out_mcf_actvol, out_mcf_fidvol, out_mcf_nmm, out_mcf_nmp, out_mcf_nem, out_mcf_nep, out_mcf_np0, out_mcf_npp, out_mcf_npm, out_mcf_nkp, out_mcf_nkm, out_mcf_nk0, out_mcf_npr, out_mcf_nne, out_mcf_nlambda, out_mcf_nsigma_p, out_mcf_nsigma_0, out_mcf_nsigma_m, out_mcf_nxi_0, out_mcf_nxi_m, out_mcf_nomega;

    out_tree->Branch("run", &out_run, "run/I");
    out_tree->Branch("sub", &out_sub, "sub/I");
    out_tree->Branch("evt", &out_evt, "evt/I");
    out_tree->Branch("selected_fiducial", &out_selected_fiducial, "selected_fiducial/I");
    out_tree->Branch("selected_topology", &out_selected_topology, "selected_topology/I");
    out_tree->Branch("reco_nu_vtx_x", &out_reco_nu_vtx_x, "reco_nu_vtx_x/F");
    out_tree->Branch("reco_nu_vtx_y", &out_reco_nu_vtx_y, "reco_nu_vtx_y/F");
    out_tree->Branch("reco_nu_vtx_z", &out_reco_nu_vtx_z, "reco_nu_vtx_z/F");
    out_tree->Branch("nhits_u", &out_nhits_u, "nhits_u/I");
    out_tree->Branch("nhit_v", &out_nhit_v, "nhit_v/I");
    out_tree->Branch("nhit_w", &out_nhit_w, "nhit_w/I");
    out_tree->Branch("charge_u", &out_charge_u, "charge_u/F");
    out_tree->Branch("charge_v", &out_charge_v, "charge_v/F");
    out_tree->Branch("charge_w", &out_charge_w, "charge_w/F");
    out_tree->Branch("wirerange_u", &out_wirerange_u, "wirerange_u/F");
    out_tree->Branch("wirerange_v", &out_wirerange_v, "wirerange_v/F");
    out_tree->Branch("wirerange_w", &out_wirerange_w, "wirerange_w/F");
    out_tree->Branch("timerange_u", &out_timerange_u, "timerange_u/F");
    out_tree->Branch("timerange_v", &out_timerange_v, "timerange_v/F");
    out_tree->Branch("timerange_w", &out_timerange_w, "timerange_w/F");
    out_tree->Branch("event_category", &out_event_category, "event_category/I");
    out_tree->Branch("mcf_nu_e", &out_mcf_nu_e, "mcf_nu_e/F");
    out_tree->Branch("mcf_lep_e", &out_mcf_lep_e, "mcf_lep_e/F");
    out_tree->Branch("mcf_actvol", &out_mcf_actvol, "mcf_actvol/I");
    out_tree->Branch("mcf_fidvol", &out_mcf_fidvol, "mcf_fidvol/I");
    out_tree->Branch("mcf_nmm", &out_mcf_nmm, "mcf_nmm/I");
    out_tree->Branch("mcf_nmp", &out_mcf_nmp, "mcf_nmp/I");
    out_tree->Branch("mcf_nem", &out_mcf_nem, "mcf_nem/I");
    out_tree->Branch("mcf_nep", &out_mcf_nep, "mcf_nep/I");
    out_tree->Branch("mcf_np0", &out_mcf_np0, "mcf_np0/I");
    out_tree->Branch("mcf_npp", &out_mcf_npp, "mcf_npp/I");
    out_tree->Branch("mcf_npm", &out_mcf_npm, "mcf_npm/I");
    out_tree->Branch("mcf_nkp", &out_mcf_nkp, "mcf_nkp/I");
    out_tree->Branch("mcf_nkm", &out_mcf_nkm, "mcf_nkm/I");
    out_tree->Branch("mcf_nk0", &out_mcf_nk0, "mcf_nk0/I");
    out_tree->Branch("mcf_npr", &out_mcf_npr, "mcf_npr/I");
    out_tree->Branch("mcf_nne", &out_mcf_nne, "mcf_nne/I");
    out_tree->Branch("mcf_nlambda", &out_mcf_nlambda, "mcf_nlambda/I");
    out_tree->Branch("mcf_nsigma_p", &out_mcf_nsigma_p, "mcf_nsigma_p/I");
    out_tree->Branch("mcf_nsigma_0", &out_mcf_nsigma_0, "mcf_nsigma_0/I");
    out_tree->Branch("mcf_nsigma_m", &out_mcf_nsigma_m, "mcf_nsigma_m/I");
    out_tree->Branch("mcf_nxi_0", &out_mcf_nxi_0, "mcf_nxi_0/I");
    out_tree->Branch("mcf_nxi_m", &out_mcf_nxi_m, "mcf_nxi_m/I");
    out_tree->Branch("mcf_nomega", &out_mcf_nomega, "mcf_nomega/I");

    for (const auto& event : unique_events) {
        int run = std::get<0>(event);
        int sub = std::get<1>(event);
        int evt = std::get<2>(event);
        out_run = run;
        out_sub = sub;
        out_evt = evt;

        Long64_t minor = static_cast<Long64_t>(sub) * multiplier + static_cast<Long64_t>(evt);
        Long64_t fid_entry = fiducial_ch.GetEntryNumberWithIndex(static_cast<Long64_t>(run), minor);
        if (fid_entry >= 0) {
            fiducial_ch.GetEntry(fid_entry);
            out_selected_fiducial = fid_selected;
            out_reco_nu_vtx_x = fid_reco_nu_vtx_x;
            out_reco_nu_vtx_y = fid_reco_nu_vtx_y;
            out_reco_nu_vtx_z = fid_reco_nu_vtx_z;
            out_event_category = fid_event_category;
            out_mcf_nu_e = fid_mcf_nu_e;
            out_mcf_lep_e = fid_mcf_lep_e;
            out_mcf_actvol = fid_mcf_actvol;
            out_mcf_fidvol = fid_mcf_fidvol;
            out_mcf_nmm = fid_mcf_nmm;
            out_mcf_nmp = fid_mcf_nmp;
            out_mcf_nem = fid_mcf_nem;
            out_mcf_nep = fid_mcf_nep;
            out_mcf_np0 = fid_mcf_np0;
            out_mcf_npp = fid_mcf_npp;
            out_mcf_npm = fid_mcf_npm;
            out_mcf_nkp = fid_mcf_nkp;
            out_mcf_nkm = fid_mcf_nkm;
            out_mcf_nk0 = fid_mcf_nk0;
            out_mcf_npr = fid_mcf_npr;
            out_mcf_nne = fid_mcf_nne;
            out_mcf_nlambda = fid_mcf_nlambda;
            out_mcf_nsigma_p = fid_mcf_nsigma_p;
            out_mcf_nsigma_0 = fid_mcf_nsigma_0;
            out_mcf_nsigma_m = fid_mcf_nsigma_m;
            out_mcf_nxi_0 = fid_mcf_nxi_0;
            out_mcf_nxi_m = fid_mcf_nxi_m;
            out_mcf_nomega = fid_mcf_nomega;
        } else {
            out_selected_fiducial = 0;
            out_reco_nu_vtx_x = -1;
            out_reco_nu_vtx_y = -1;
            out_reco_nu_vtx_z = -1;
        }

        Long64_t top_entry = topology_ch.GetEntryNumberWithIndex(static_cast<Long64_t>(run), minor);
        if (top_entry >= 0) {
            topology_ch.GetEntry(top_entry);
            out_selected_topology = top_selected;
            out_nhits_u = top_nhits_u;
            out_nhit_v = top_nhit_v;
            out_nhit_w = top_nhit_w;
            out_charge_u = top_charge_u;
            out_charge_v = top_charge_v;
            out_charge_w = top_charge_w;
            out_wirerange_u = top_wirerange_u;
            out_wirerange_v = top_wirerange_v;
            out_wirerange_w = top_wirerange_w;
            out_timerange_u = top_timerange_u;
            out_timerange_v = top_timerange_v;
            out_timerange_w = top_timerange_w;
            if (fid_entry < 0) {
                out_event_category = top_event_category;
                out_mcf_nu_e = top_mcf_nu_e;
                out_mcf_lep_e = top_mcf_lep_e;
                out_mcf_actvol = top_mcf_actvol;
                out_mcf_fidvol = top_mcf_fidvol;
                out_mcf_nmm = top_mcf_nmm;
                out_mcf_nmp = top_mcf_nmp;
                out_mcf_nem = top_mcf_nem;
                out_mcf_nep = top_mcf_nep;
                out_mcf_np0 = top_mcf_np0;
                out_mcf_npp = top_mcf_npp;
                out_mcf_npm = top_mcf_npm;
                out_mcf_nkp = top_mcf_nkp;
                out_mcf_nkm = top_mcf_nkm;
                out_mcf_nk0 = top_mcf_nk0;
                out_mcf_npr = top_mcf_npr;
                out_mcf_nne = top_mcf_nne;
                out_mcf_nlambda = top_mcf_nlambda;
                out_mcf_nsigma_p = top_mcf_nsigma_p;
                out_mcf_nsigma_0 = top_mcf_nsigma_0;
                out_mcf_nsigma_m = top_mcf_nsigma_m;
                out_mcf_nxi_0 = top_mcf_nxi_0;
                out_mcf_nxi_m = top_mcf_nxi_m;
                out_mcf_nomega = top_mcf_nomega;
            }
        } else {
            out_selected_topology = 0;
            out_nhits_u = -1;
            out_nhit_v = -1;
            out_nhit_w = -1;
            out_charge_u = -1;
            out_charge_v = -1;
            out_charge_w = -1;
            out_wirerange_u = -1;
            out_wirerange_v = -1;
            out_wirerange_w = -1;
            out_timerange_u = -1;
            out_timerange_v = -1;
            out_timerange_w = -1;
        }

        out_tree->Fill();
    }

    TChain subrun_ch("topologyTrainingFilter/SubRun");
    for (const auto& f_name : in_file_name) {
        subrun_ch.Add(f_name.c_str());
    }

    float pot;
    float summed_pot = 0;
    subrun_ch.SetBranchAddress("pot", &pot);
    for (Long64_t i = 0; i < subrun_ch.GetEntries(); ++i) {
        subrun_ch.GetEntry(i);
        summed_pot += pot;
    }

    TParameter<float>* summed_pot_param = new TParameter<float>("summed_pot", summed_pot);
    summed_pot_param->Write();

    out_tree->Write();
    out_file->Close();
    delete out_file;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "-- Usage: analyse INPUT_FILE OUTPUT_FILE\n";
        return 1;
    }

    std::string input_file_name(argv[1]);
    std::string output_file_name(argv[2]);

    std::vector<std::string> in_files = {input_file_name};
    analyse(in_files, output_file_name);

    return 0;
}