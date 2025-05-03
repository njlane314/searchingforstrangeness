#include <iostream>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "EventCategory.hh"
#include "MyPointer.hh"
#include "TreeUtils.hh"

class AnalysisEvent {
public:
    AnalysisEvent() {}
    ~AnalysisEvent() {}

    void set_branch_addresses(TTree& tree);
    void set_output_branch_addresses(TTree& out_tree, bool create = false);
    void categorise_event();
    void passes_fiducial_selection();

    bool _is_mc;
    float _mcf_nu_e = 0;
    float _mcf_lep_e = 0;
    int _mcf_nu_pdg = 0;
    int _mcf_ccnc = 0;
    int _mcf_mode = 0;
    int _mcf_interaction = 0;
    int _mcf_target_nucleus = 0;
    int _mcf_hit_nucleon = 0;
    float _mcf_W = 0;
    float _mcf_X = 0;
    float _mcf_Y = 0;
    float _mcf_QSqr = 0;
    float _mcf_px = 0;
    float _mcf_py = 0;
    float _mcf_pz = 0;
    float _mcf_theta = 0;
    float _mcf_pt = 0;
    int _mcf_actvol = 0;
    int _mcf_fidvol = 0;
    int _mcf_nmm = 0;
    int _mcf_nmp = 0;
    int _mcf_nem = 0;
    int _mcf_nep = 0;
    int _mcf_np0 = 0;
    int _mcf_npp = 0;
    int _mcf_npm = 0;
    int _mcf_nkp = 0;
    int _mcf_nkm = 0;
    int _mcf_nk0 = 0;
    int _mcf_npr = 0;
    int _mcf_nne = 0;
    int _mcf_nlambda = 0;
    int _mcf_nsigma_p = 0;
    int _mcf_nsigma_0 = 0;
    int _mcf_nsigma_m = 0;
    int _mcf_nxi_0 = 0;
    int _mcf_nxi_m = 0;
    int _mcf_nomega = 0;

    bool in_fiducial; 
    int m_nhit_u = 0;
    int m_nhits_v = 0;
    int m_nhits_w = 0;
    float m_charge_u = 0;
    float m_charge_v = 0;
    float m_charge_w = 0;
    float m_wirerange_u = 0;
    float m_wirerange_v = 0;
    float m_wirerange_w = 0;
    float m_timerange_u = 0;
    float m_timerange_v = 0;
    float m_timerange_w = 0;
 
    MyPointer<std::vector<std::vector<float>>> calo_pixels_u;
    MyPointer<std::vector<std::vector<float>>> calo_pixels_v;
    MyPointer<std::vector<std::vector<float>>> calo_pixels_w;
    MyPointer<std::vector<std::vector<float>>> reco_pixels_u;
    MyPointer<std::vector<std::vector<float>>> reco_pixels_v;
    MyPointer<std::vector<std::vector<float>>> reco_pixels_w;
    MyPointer<std::vector<std::vector<float>>> label_pixels_u;
    MyPointer<std::vector<std::vector<float>>> label_pixels_v;
    MyPointer<std::vector<std::vector<float>>> label_pixels_w;   

    int _event_category = static_cast<int>(EventCategory::k_unknown);

    bool sel_reco_vertex_in_fv;  
};

void AnalysisEvent::set_branch_addresses(TTree& tree) {
    tree.SetBranchAddress("is_mc", &_is_mc);
    tree.SetBranchAddress("mcf_nu_e", &_mcf_nu_e);
    tree.SetBranchAddress("mcf_lep_e", &_mcf_lep_e);
    tree.SetBranchAddress("mcf_nu_pdg", &_mcf_nu_pdg);
    tree.SetBranchAddress("mcf_ccnc", &_mcf_ccnc);
    tree.SetBranchAddress("mcf_mode", &_mcf_mode);
    tree.SetBranchAddress("mcf_interaction", &_mcf_interaction);
    tree.SetBranchAddress("mcf_target_nucleus", &_mcf_target_nucleus);
    tree.SetBranchAddress("mcf_hit_nucleon", &_mcf_hit_nucleon);
    tree.SetBranchAddress("mcf_W", &_mcf_W);
    tree.SetBranchAddress("mcf_X", &_mcf_X);
    tree.SetBranchAddress("mcf_Y", &_mcf_Y);
    tree.SetBranchAddress("mcf_QSqr", &_mcf_QSqr);
    tree.SetBranchAddress("mcf_px", &_mcf_px);
    tree.SetBranchAddress("mcf_py", &_mcf_py);
    tree.SetBranchAddress("mcf_pz", &_mcf_pz);
    tree.SetBranchAddress("mcf_theta", &_mcf_theta);
    tree.SetBranchAddress("mcf_pt", &_mcf_pt);
    tree.SetBranchAddress("mcf_actvol", &_mcf_actvol);
    tree.SetBranchAddress("mcf_fidvol", &_mcf_fidvol);
    tree.SetBranchAddress("mcf_nmm", &_mcf_nmm);
    tree.SetBranchAddress("mcf_nmp", &_mcf_nmp);
    tree.SetBranchAddress("mcf_nem", &_mcf_nem);
    tree.SetBranchAddress("mcf_nep", &_mcf_nep);
    tree.SetBranchAddress("mcf_np0", &_mcf_np0);
    tree.SetBranchAddress("mcf_npp", &_mcf_npp);
    tree.SetBranchAddress("mcf_npm", &_mcf_npm);
    tree.SetBranchAddress("mcf_nkp", &_mcf_nkp);
    tree.SetBranchAddress("mcf_nkm", &_mcf_nkm);
    tree.SetBranchAddress("mcf_nk0", &_mcf_nk0);
    tree.SetBranchAddress("mcf_npr", &_mcf_npr);
    tree.SetBranchAddress("mcf_nne", &_mcf_nne);
    tree.SetBranchAddress("mcf_nlambda", &_mcf_nlambda);
    tree.SetBranchAddress("mcf_nsigma_p", &_mcf_nsigma_p);
    tree.SetBranchAddress("mcf_nsigma_0", &_mcf_nsigma_0);
    tree.SetBranchAddress("mcf_nsigma_m", &_mcf_nsigma_m);
    tree.SetBranchAddress("mcf_nxi_0", &_mcf_nxi_0);
    tree.SetBranchAddress("mcf_nxi_m", &_mcf_nxi_m);
    tree.SetBranchAddress("mcf_nomega", &_mcf_nomega);

    tree.SetBranchAddress("in_fiducial", &in_fiducial);
    tree.SetBranchAddress("nhit_u", &m_nhit_u);
    tree.SetBranchAddress("nhit_v", &m_nhits_v);
    tree.SetBranchAddress("nhit_w", &m_nhits_w);
    tree.SetBranchAddress("charge_u", &m_charge_u);
    tree.SetBranchAddress("charge_v", &m_charge_v);
    tree.SetBranchAddress("charge_w", &m_charge_w);
    tree.SetBranchAddress("wirerange_u", &m_wirerange_u);
    tree.SetBranchAddress("wirerange_v", &m_wirerange_v);
    tree.SetBranchAddress("wirerange_w", &m_wirerange_w);
    tree.SetBranchAddress("timerange_u", &m_timerange_u);
    tree.SetBranchAddress("timerange_v", &m_timerange_v);
    tree.SetBranchAddress("timerange_w", &m_timerange_w);

    set_object_input_branch_address(tree, "calo_pixels_u", calo_pixels_u);
    set_object_input_branch_address(tree, "calo_pixels_v", calo_pixels_v);
    set_object_input_branch_address(tree, "calo_pixels_w", calo_pixels_w);
    set_object_input_branch_address(tree, "reco_pixels_u", reco_pixels_u);
    set_object_input_branch_address(tree, "reco_pixels_v", reco_pixels_v);
    set_object_input_branch_address(tree, "reco_pixels_w", reco_pixels_w);
    set_object_input_branch_address(tree, "label_pixels_u", label_pixels_u);
    set_object_input_branch_address(tree, "label_pixels_v", label_pixels_v);
    set_object_input_branch_address(tree, "label_pixels_w", label_pixels_w);
}

void AnalysisEvent::set_output_branch_addresses(TTree& out_tree, bool create) {
    set_output_branch_address(out_tree, "is_mc", &_is_mc, create, "is_mc/O");
    set_output_branch_address(out_tree, "mcf_nu_e", &_mcf_nu_e, create, "mcf_nu_e/F");
    set_output_branch_address(out_tree, "mcf_lep_e", &_mcf_lep_e, create, "mcf_lep_e/F");
    set_output_branch_address(out_tree, "mcf_nu_pdg", &_mcf_nu_pdg, create, "mcf_nu_pdg/I");
    set_output_branch_address(out_tree, "mcf_ccnc", &_mcf_ccnc, create, "mcf_ccnc/I");
    set_output_branch_address(out_tree, "mcf_mode", &_mcf_mode, create, "mcf_mode/I");
    set_output_branch_address(out_tree, "mcf_interaction", &_mcf_interaction, create, "mcf_interaction/I");
    set_output_branch_address(out_tree, "mcf_target_nucleus", &_mcf_target_nucleus, create, "mcf_target_nucleus/I");
    set_output_branch_address(out_tree, "mcf_hit_nucleon", &_mcf_hit_nucleon, create, "mcf_hit_nucleon/I");
    set_output_branch_address(out_tree, "mcf_W", &_mcf_W, create, "mcf_W/F");
    set_output_branch_address(out_tree, "mcf_X", &_mcf_X, create, "mcf_X/F");
    set_output_branch_address(out_tree, "mcf_Y", &_mcf_Y, create, "mcf_Y/F");
    set_output_branch_address(out_tree, "mcf_QSqr", &_mcf_QSqr, create, "mcf_QSqr/F");
    set_output_branch_address(out_tree, "mcf_px", &_mcf_px, create, "mcf_px/F");
    set_output_branch_address(out_tree, "mcf_py", &_mcf_py, create, "mcf_py/F");
    set_output_branch_address(out_tree, "mcf_pz", &_mcf_pz, create, "mcf_pz/F");
    set_output_branch_address(out_tree, "mcf_theta", &_mcf_theta, create, "mcf_theta/F");
    set_output_branch_address(out_tree, "mcf_pt", &_mcf_pt, create, "mcf_pt/F");
    set_output_branch_address(out_tree, "mcf_actvol", &_mcf_actvol, create, "mcf_actvol/I");
    set_output_branch_address(out_tree, "mcf_fidvol", &_mcf_fidvol, create, "mcf_fidvol/I");
    set_output_branch_address(out_tree, "mcf_nmm", &_mcf_nmm, create, "mcf_nmm/I");
    set_output_branch_address(out_tree, "mcf_nmp", &_mcf_nmp, create, "mcf_nmp/I");
    set_output_branch_address(out_tree, "mcf_nem", &_mcf_nem, create, "mcf_nem/I");
    set_output_branch_address(out_tree, "mcf_nep", &_mcf_nep, create, "mcf_nep/I");
    set_output_branch_address(out_tree, "mcf_np0", &_mcf_np0, create, "mcf_np0/I");
    set_output_branch_address(out_tree, "mcf_npp", &_mcf_npp, create, "mcf_npp/I");
    set_output_branch_address(out_tree, "mcf_npm", &_mcf_npm, create, "mcf_npm/I");
    set_output_branch_address(out_tree, "mcf_nkp", &_mcf_nkp, create, "mcf_nkp/I");
    set_output_branch_address(out_tree, "mcf_nkm", &_mcf_nkm, create, "mcf_nkm/I");
    set_output_branch_address(out_tree, "mcf_nk0", &_mcf_nk0, create, "mcf_nk0/I");
    set_output_branch_address(out_tree, "mcf_npr", &_mcf_npr, create, "mcf_npr/I");
    set_output_branch_address(out_tree, "mcf_nne", &_mcf_nne, create, "mcf_nne/I");
    set_output_branch_address(out_tree, "mcf_nlambda", &_mcf_nlambda, create, "mcf_nlambda/I");
    set_output_branch_address(out_tree, "mcf_nsigma_p", &_mcf_nsigma_p, create, "mcf_nsigma_p/I");
    set_output_branch_address(out_tree, "mcf_nsigma_0", &_mcf_nsigma_0, create, "mcf_nsigma_0/I");
    set_output_branch_address(out_tree, "mcf_nsigma_m", &_mcf_nsigma_m, create, "mcf_nsigma_m/I");
    set_output_branch_address(out_tree, "mcf_nxi_0", &_mcf_nxi_0, create, "mcf_nxi_0/I");
    set_output_branch_address(out_tree, "mcf_nxi_m", &_mcf_nxi_m, create, "mcf_nxi_m/I");
    set_output_branch_address(out_tree, "mcf_nomega", &_mcf_nomega, create, "mcf_nomega/I");

    set_output_branch_address(out_tree, "in_fiducial", &in_fiducial, create, "in_fiducial/O");
    set_output_branch_address(out_tree, "nhit_u", &m_nhit_u, create, "nhit_u/I");
    set_output_branch_address(out_tree, "nhit_v", &m_nhits_v, create, "nhit_v/I");
    set_output_branch_address(out_tree, "nhit_w", &m_nhits_w, create, "nhit_w/I");
    set_output_branch_address(out_tree, "charge_u", &m_charge_u, create, "charge_u/F");
    set_output_branch_address(out_tree, "charge_v", &m_charge_v, create, "charge_v/F");
    set_output_branch_address(out_tree, "charge_w", &m_charge_w, create, "charge_w/F");
    set_output_branch_address(out_tree, "wirerange_u", &m_wirerange_u, create, "wirerange_u/F");
    set_output_branch_address(out_tree, "wirerange_v", &m_wirerange_v, create, "wirerange_v/F");
    set_output_branch_address(out_tree, "wirerange_w", &m_wirerange_w, create, "wirerange_w/F");
    set_output_branch_address(out_tree, "timerange_u", &m_timerange_u, create, "timerange_u/F");
    set_output_branch_address(out_tree, "timerange_v", &m_timerange_v, create, "timerange_v/F");
    set_output_branch_address(out_tree, "timerange_w", &m_timerange_w, create, "timerange_w/F");
   
    set_object_output_branch_address<std::vector<std::vector<float>>>(out_tree, "calo_pixels_u", calo_pixels_u, create);
    set_object_output_branch_address<std::vector<std::vector<float>>>(out_tree, "calo_pixels_v", calo_pixels_v, create);
    set_object_output_branch_address<std::vector<std::vector<float>>>(out_tree, "calo_pixels_w", calo_pixels_w, create);
    set_object_output_branch_address<std::vector<std::vector<float>>>(out_tree, "reco_pixels_u", reco_pixels_u, create);
    set_object_output_branch_address<std::vector<std::vector<float>>>(out_tree, "reco_pixels_v", reco_pixels_v, create);
    set_object_output_branch_address<std::vector<std::vector<float>>>(out_tree, "reco_pixels_w", reco_pixels_w, create);
    set_object_output_branch_address<std::vector<std::vector<float>>>(out_tree, "label_pixels_u", label_pixels_u, create);
    set_object_output_branch_address<std::vector<std::vector<float>>>(out_tree, "label_pixels_v", label_pixels_v, create);
    set_object_output_branch_address<std::vector<std::vector<float>>>(out_tree, "label_pixels_w", label_pixels_w, create);

    set_output_branch_address(out_tree, "event_category", &_event_category, create, "event_category/I");

    set_output_branch_address(out_tree, "sel_reco_vertex_in_fv", &sel_reco_vertex_in_fv, create, "sel_reco_vertex_in_fv/O");
}

void AnalysisEvent::categorise_event() {
    if (!_is_mc || !_mcf_actvol) {
        _event_category = static_cast<int>(EventCategory::k_unknown);
        return;
    }

    if (_mcf_ccnc == 1) {
        _event_category = static_cast<int>(EventCategory::k_nc);
    } else {
        int lepton_pdg = std::abs(_mcf_nu_pdg);
        if (lepton_pdg == 11) {
            _event_category = static_cast<int>(EventCategory::k_nu_e_cc);
        } else if (lepton_pdg == 13) {
            bool has_strange_hadron = (_mcf_nkp > 0 || _mcf_nkm > 0 || _mcf_nk0 > 0 ||
                                       _mcf_nlambda > 0 || _mcf_nsigma_p > 0 || _mcf_nsigma_0 > 0 ||
                                       _mcf_nsigma_m > 0 || _mcf_nxi_0 > 0 || _mcf_nxi_m > 0 ||
                                       _mcf_nomega > 0);
            bool has_protons = (_mcf_npr > 0);
            bool has_pions = (_mcf_npp > 0 || _mcf_npm > 0);
            bool has_pi0 = (_mcf_np0 > 0);
            if (has_strange_hadron) {
                _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_with_strange);
            } else if (has_protons && !has_pions && !has_pi0) {
                _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_with_protons);
            } else if (!has_protons && has_pions && !has_pi0) {
                _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_with_pions);
            } else if (has_protons && has_pions && !has_pi0) {
                _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_with_protons_pions);
            } else {
                _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_other);
            }
        } else {
            _event_category = static_cast<int>(EventCategory::k_other);
        }
    }
}

void AnalysisEvent::passes_fiducial_selection() {
    sel_reco_vertex_in_fv = in_fiducial;
}

void analyse(const std::vector<std::string>& in_file_names, const std::string& output_filename) {
    TChain events_ch("strangenessFilter/EventSelectionFilter");
    TChain subruns_ch("strangenessFilter/SubRun");

    for (const auto& f_name : in_file_names) {
        events_ch.Add(f_name.c_str());
        subruns_ch.Add(f_name.c_str());
    }

    TFile* out_file = new TFile(output_filename.c_str(), "recreate");
    out_file->cd();
    TTree* out_tree = new TTree("analysis_tree", "");

    float pot;
    float summed_pot = 0.;
    bool has_pot_branch = (subruns_ch.GetBranch("pot") != nullptr);
    if (has_pot_branch) {
        subruns_ch.SetBranchAddress("pot", &pot);
        for (int se = 0; se < subruns_ch.GetEntries(); ++se) {
            subruns_ch.GetEntry(se);
            summed_pot += pot;
        }
    }
    TParameter<float>* summed_pot_param = new TParameter<float>("summed_pot", summed_pot);
    summed_pot_param->Write();

    bool created_output_branches = false;
    long events_entry = 0;
    while (true) {
        if (events_entry % 1000 == 0) {
            std::cout << "Processing event #" << events_entry << '\n';
        }

        AnalysisEvent cur_event;
        cur_event.set_branch_addresses(events_ch);

        int local_entry = events_ch.LoadTree(events_entry);
        if (local_entry < 0) break;

        events_ch.GetEntry(events_entry);

        if (!created_output_branches) {
            cur_event.set_output_branch_addresses(*out_tree, true);
            created_output_branches = true;
        } else {
            cur_event.set_output_branch_addresses(*out_tree);
        }

        cur_event.categorise_event();
        cur_event.passes_fiducial_selection();
        out_tree->Fill();
        ++events_entry;
    }

    out_tree->Write();
    out_file->Close();
    delete out_file;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: analyzer INPUT_FILE OUTPUT_FILE\n";
        return 1;
    }

    std::string input_file_name(argv[1]);
    std::string output_file_name(argv[2]);

    std::vector<std::string> in_files = {input_file_name};
    analyse(in_files, output_file_name);

    return 0;
}