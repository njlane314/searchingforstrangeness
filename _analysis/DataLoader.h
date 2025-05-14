#pragma once

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <functional>

#include "TFile.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/VectorUtil.h"

namespace AnalysisFramework {

enum class NTupleFileType {
    kUnknown = 0, kDataBNB, kDataNuMI, kBackgroundEXT_BNB, kBackgroundEXT_NuMI,
    kMCInclusiveBNB, kMCInclusiveNuMI, kMCStrangenessSignalBNB, kMCStrangenessSignalNuMI,
    kMCDirtBNB, kMCDirtNuMI, kMCGeneric, kDetVarCV, kDetVarLYAttenuation, kDetVarLYDown,
    kDetVarLYRayleigh, kDetVarRecomb2, kDetVarSCE, kDetVarWireModAngleXZ, kDetVarWireModAngleYZ,
};

inline bool IsBeamOnData(NTupleFileType type) { return type == NTupleFileType::kDataBNB || type == NTupleFileType::kDataNuMI; }
inline bool IsBeamOffBackground(NTupleFileType type) { return type == NTupleFileType::kBackgroundEXT_BNB || type == NTupleFileType::kBackgroundEXT_NuMI; }
inline bool IsMonteCarlo(NTupleFileType type) { return !IsBeamOnData(type) && !IsBeamOffBackground(type) && type != NTupleFileType::kUnknown; }
inline bool IsIntrinsicStrangenessSignalMC(NTupleFileType type) { return type == NTupleFileType::kMCStrangenessSignalBNB || type == NTupleFileType::kMCStrangenessSignalNuMI; }
inline bool IsDetectorVariation(NTupleFileType type) {
    return type == NTupleFileType::kDetVarCV || type == NTupleFileType::kDetVarLYAttenuation ||
           type == NTupleFileType::kDetVarLYDown || type == NTupleFileType::kDetVarLYRayleigh ||
           type == NTupleFileType::kDetVarRecomb2 || type == NTupleFileType::kDetVarSCE ||
           type == NTupleFileType::kDetVarWireModAngleXZ || type == NTupleFileType::kDetVarWireModAngleYZ;
}

struct FilePathConfig {
    std::string ntuple_path_base;
    std::string analysis_config_cpp_path;
    FilePathConfig(std::string base = "/pnfs/uboone/persistent/users/YourUser/ntuples/", std::string cfg = "./config/")
        : ntuple_path_base(std::move(base)), analysis_config_cpp_path(std::move(cfg)) {}
};

struct DatasetInfo {
    std::string file_prefix; double pot = 0.0; int triggers = 0;
    std::string filter_query; NTupleFileType type = NTupleFileType::kUnknown; std::string name_in_code;
};
struct RunConfig {
    std::string run_id_key; std::string beam_type; std::string relative_path;
    std::map<std::string, DatasetInfo> datasets;
};
struct DetVarDatasetInfo {
    std::string file_prefix_variant_suffix; double pot_scaling_factor_if_different = 1.0;
    NTupleFileType base_sample_type = NTupleFileType::kUnknown; std::string name_in_code_cv_equivalent;
};
struct DetVarRunConfig {
    std::string run_id_key; std::string beam_type; std::string relative_path_detvar_base;
    std::map<std::string, std::map<std::string, DetVarDatasetInfo>> variations;
};

class VariableManager {
public:
    VariableManager() {
        base_event_id_vars_ = {"run", "sub", "evt"};
        reco_event_vars_ = {
            "category", "reco_nu_vtx_x", "reco_nu_vtx_y", "reco_nu_vtx_z", "reco_nu_vtx_sce_x", "reco_nu_vtx_sce_y", "reco_nu_vtx_sce_z",
            "slpdg", "nslice", "slnhits", "slclustfrac", "n_showers_contained", "n_tracks_contained", "nelec", "nproton", "npion",
            "nneutron", "muon_e", "electron_e", "proton_e", "pion_e", "endmuonmichel", "CosmicIP", "CosmicDirAll3D", "CosmicIPAll3D",
            "_opfilter_pe_beam", "_opfilter_pe_veto", "topological_score", "isVtxInFiducial", "theta"
        };
        mc_truth_event_vars_ = {
            "nu_pdg", "nu_ccnc", "nu_mode", "nu_interaction", "nu_e", "nu_theta", "nu_pt", "nu_target_nucleus", "nu_hit_nucleon",
            "nu_W", "nu_X", "nu_Y", "nu_QSqr", "nu_px", "nu_py", "nu_pz", "nu_vtx_x", "nu_vtx_y", "nu_vtx_z", "nu_in_fidvol",
            "mcf_nmm", "mcf_nmp", "mcf_nem", "mcf_nep", "mcf_np0", "mcf_npp", "mcf_npm", "mcf_npr", "mcf_nne",
            "mcf_nkp", "mcf_nkm", "mcf_nk0", "mcf_nlambda", "mcf_nsigma_p", "mcf_nsigma_0", "mcf_nsigma_m",
            "mcf_nxi_0", "mcf_nxi_m", "mcf_nomega"
        };
        reco_track_vars_ = {
            "trk_pfp_id_v", "trk_start_x_v", "trk_start_y_v", "trk_start_z_v", "trk_sce_start_x_v", "trk_sce_start_y_v",
            "trk_sce_start_z_v", "trk_end_x_v", "trk_end_y_v", "trk_end_z_v", "trk_sce_end_x_v", "trk_sce_end_y_v",
            "trk_sce_end_z_v", "trk_dir_x_v", "trk_dir_y_v", "trk_dir_z_v", "trk_len_v", "trk_distance_v", "trk_theta_v",
            "trk_phi_v", "trk_bragg_p_v", "trk_bragg_mu_v", "trk_bragg_pion_v", "trk_bragg_mip_v", "trk_pid_chipr_v",
            "trk_pid_chika_v", "trk_pid_chipi_v", "trk_pid_chimu_v", "trk_pida_v", "trk_llr_pid_u_v", "trk_llr_pid_v_v",
            "trk_llr_pid_y_v", "trk_llr_pid_v", "trk_llr_pid_score_v", "trk_mcs_muon_mom_v", "trk_range_muon_mom_v",
            "trk_energy_proton_v", "trk_energy_muon_v", "trk_calo_energy_u_v", "trk_calo_energy_v_v", "trk_calo_energy_y_v",
            "trk_trunk_dEdx_u_v", "trk_trunk_dEdx_v_v", "trk_trunk_dEdx_y_v", "trk_trunk_rr_dEdx_u_v", "trk_trunk_rr_dEdx_v_v",
            "trk_trunk_rr_dEdx_y_v", "trk_nhits_u_v", "trk_nhits_v_v", "trk_nhits_y_v", "trk_avg_deflection_mean_v",
            "trk_avg_deflection_stdev_v", "trk_avg_deflection_separation_mean_v", "trk_end_spacepoints_v", "trk_bt_pdg_v",
            "trk_score_v", "pfp_generation_v"
        };
        nominal_mc_weights_ = {"weightSpline", "weightTune", "weightSplineTimesTune", "ppfx_cv", "nu_decay_mode"};
        systematic_knob_weights_ = {
            "knobRPAup", "knobRPAdn", "knobCCMECup", "knobCCMECdn", "knobAxFFCCQEup", "knobAxFFCCQEdn",
            "knobVecFFCCQEup", "knobVecFFCCQEdn", "knobDecayAngMECup", "knobDecayAngMECdn",
            "knobThetaDelta2Npiup", "knobThetaDelta2Npidn", "knobThetaDelta2NRadup", "knobThetaDelta2NRaddn",
            "knobNormCCCOHup", "knobNormCCCOHdn", "knobNormNCCOHup", "knobNormNCCOHdn",
            "knobxsr_scc_Fv3up", "knobxsr_scc_Fv3dn", "knobxsr_scc_Fa3up", "knobxsr_scc_Fa3dn", "RootinoFix"
        };
        multi_universe_weights_ = {"weightsGenie", "weightsFlux", "weightsReint", "weightsPPFX"};
        intrinsic_strangeness_specific_event_weight_branch_ = {"is_channel_correction_weight"}; // Example branch name
    }

    std::vector<std::string> GetVariables(
        bool load_reco_event_info, bool load_reco_track_info, bool load_mc_truth_info,
        bool load_mc_weights_and_systematics, bool load_specific_is_event_weight_branch,
        NTupleFileType sample_type
    ) const {
        std::set<std::string> vars_set;
        vars_set.insert(base_event_id_vars_.begin(), base_event_id_vars_.end());
        if (load_reco_event_info) vars_set.insert(reco_event_vars_.begin(), reco_event_vars_.end());
        if (load_reco_track_info) vars_set.insert(reco_track_vars_.begin(), reco_track_vars_.end());
        if (load_mc_truth_info && IsMonteCarlo(sample_type)) {
            vars_set.insert(mc_truth_event_vars_.begin(), mc_truth_event_vars_.end());
        }
        if (load_mc_weights_and_systematics && IsMonteCarlo(sample_type)) {
            vars_set.insert(nominal_mc_weights_.begin(), nominal_mc_weights_.end());
            vars_set.insert(systematic_knob_weights_.begin(), systematic_knob_weights_.end());
            vars_set.insert(multi_universe_weights_.begin(), multi_universe_weights_.end());
        }
        if (load_specific_is_event_weight_branch && IsIntrinsicStrangenessSignalMC(sample_type)) {
            vars_set.insert(intrinsic_strangeness_specific_event_weight_branch_.begin(), intrinsic_strangeness_specific_event_weight_branch_.end());
        }
        return std::vector<std::string>(vars_set.begin(), vars_set.end());
    }
private:
    std::vector<std::string> base_event_id_vars_, reco_event_vars_, mc_truth_event_vars_, reco_track_vars_;
    std::vector<std::string> nominal_mc_weights_, systematic_knob_weights_, multi_universe_weights_;
    std::vector<std::string> intrinsic_strangeness_specific_event_weight_branch_;
};

class ConfigurationManager {
public:
    ConfigurationManager(FilePathConfig file_paths) : file_paths_(std::move(file_paths)) {}
    void PopulateConfigurations(const std::string& config_set_name = "default"); // To be defined in user's .cpp
    const RunConfig& GetRunConfig(const std::string& run_id_key) const {
        auto it = run_configs_.find(run_id_key);
        if (it == run_configs_.end()) throw std::runtime_error("RunConfig not found: " + run_id_key);
        return it->second;
    }
    const DetVarRunConfig& GetDetVarRunConfig(const std::string& run_id_key) const {
        auto it = det_var_run_configs_.find(run_id_key);
        if (it == det_var_run_configs_.end()) throw std::runtime_error("DetVarRunConfig not found: " + run_id_key);
        return it->second;
    }
    const FilePathConfig& GetFilePaths() const { return file_paths_; }
    void AddRunConfig(const RunConfig& config) { run_configs_[config.run_id_key] = config; }
    void AddDetVarRunConfig(const DetVarRunConfig& config) { det_var_run_configs_[config.run_id_key] = config; }
private:
    FilePathConfig file_paths_;
    std::map<std::string, RunConfig> run_configs_;
    std::map<std::string, DetVarRunConfig> det_var_run_configs_;
};

namespace DataLoader {

template<typename T_vec, typename T_val = typename T_vec::value_type>
T_val GetElementFromVector(const T_vec& vec, int index, T_val default_val = T_val{}) {
    if (index >= 0 && static_cast<size_t>(index) < vec.size()) return vec[index];
    return default_val;
}
int GetIndexFromVectorSort(const ROOT::RVec<float>& values_vec, const ROOT::RVec<bool>& mask_vec, int n_th_idx = 0, bool asc = false) {
    if (values_vec.empty() || (!mask_vec.empty() && values_vec.size()!=mask_vec.size())) return -1;
    std::vector<std::pair<float,int>> masked;
    if(!mask_vec.empty()) { for(size_t i=0;i<values_vec.size();++i) if(mask_vec[i]) masked.push_back({values_vec[i],(int)i}); }
    else { for(size_t i=0;i<values_vec.size();++i) masked.push_back({values_vec[i],(int)i}); }
    if (masked.empty()||n_th_idx<0||(size_t)n_th_idx>=masked.size()) return -1;
    auto comp = asc ? [](const auto&a,const auto&b){return a.first<b.first;} : [](const auto&a,const auto&b){return a.first>b.first;};
    std::nth_element(masked.begin(),masked.begin()+n_th_idx,masked.end(),comp);
    return masked[n_th_idx].second;
}

ROOT::RDF::RInterface AddEventCategories(ROOT::RDF::RInterface df, const std::string& dataset_key_for_event_cat_logic) {
    auto df_with_defs = df;
    bool is_mc_for_cat = (dataset_key_for_event_cat_logic != "data" && dataset_key_for_event_cat_logic != "ext_background");

    if (is_mc_for_cat) {
        bool has_npp=false,has_npm=false,has_npr=false;
        for(const auto&cn:df.GetColumnNames()){if(cn=="mcf_npp")has_npp=true;if(cn=="mcf_npm")has_npm=true;if(cn=="mcf_npr")has_npr=true;}
        if(has_npp&&has_npm)df_with_defs=df_with_defs.Define("mc_n_charged_pions_true","mcf_npp+mcf_npm");else df_with_defs=df_with_defs.Define("mc_n_charged_pions_true","return -1;");
        if(has_npr)df_with_defs=df_with_defs.Define("mc_n_protons_true","mcf_npr");else df_with_defs=df_with_defs.Define("mc_n_protons_true","return -1;");
    } else {
        df_with_defs=df_with_defs.Define("mc_n_charged_pions_true","return -1;").Define("mc_n_protons_true","return -1;");
    }

    auto df_total_strange = df_with_defs.Define("total_true_strange_particles",
        [](int nkp,int nkm,int nk0,int nl,int nsp,int ns0,int nsm,int nx0,int nxm,int nom){return nkp+nkm+nk0+nl+nsp+ns0+nsm+nx0+nxm+nom;},
        {"mcf_nkp","mcf_nkm","mcf_nk0","mcf_nlambda","mcf_nsigma_p","mcf_nsigma_0","mcf_nsigma_m","mcf_nxi_0","mcf_nxi_m","mcf_nomega"});
    auto df_strange_multiplicity = df_total_strange.Define("inclusive_strangeness_multiplicity_type",
        [](int ts){if(ts==0)return 0;if(ts==1)return 1;return 2;},{"total_true_strange_particles"});

    auto df_with_event_category = df_strange_multiplicity.Define("event_category_val",
        [dataset_key_for_event_cat_logic](int nu_pdg,int ccnc,int npi_char_true,int npr_true,int str_mult){
            int cat=9999;
            if(dataset_key_for_event_cat_logic=="data")cat=0;
            else if(dataset_key_for_event_cat_logic=="ext_background")cat=1;
            else if(dataset_key_for_event_cat_logic=="dirt_mc")cat=2;
            else if(dataset_key_for_event_cat_logic == "intrinsic_strangeness_signal" || dataset_key_for_event_cat_logic == "inclusive_generic_mc"){
                bool isnumu=std::abs(nu_pdg)==14; bool isnue=std::abs(nu_pdg)==12; bool iscc=(ccnc==0); bool isnc=(ccnc==1);
                if(str_mult>0){
                    if(isnumu&&iscc)cat=(str_mult==1)?10:11; else if(isnue&&iscc)cat=(str_mult==1)?12:13;
                    else if(isnc)cat=(str_mult==1)?14:15; else cat=19;
                }else{
                    if(isnumu){
                        if(iscc){
                            if(npi_char_true==0){if(npr_true==0)cat=100;else if(npr_true==1)cat=101;else cat=102;}
                            else if(npi_char_true==1){if(npr_true==0)cat=103;else if(npr_true==1)cat=104;else cat=105;}
                            else cat=106;
                        }else{
                            if(npi_char_true==0){if(npr_true==0)cat=110;else if(npr_true==1)cat=111;else cat=112;}
                            else if(npi_char_true==1){if(npr_true==0)cat=113;else if(npr_true==1)cat=114;else cat=115;}
                            else cat=116;
                        }
                    }else if(isnue){
                        if(iscc){
                            if(npi_char_true==0){if(npr_true==0)cat=200;else if(npr_true==1)cat=201;else cat=202;}
                            else if(npi_char_true==1){if(npr_true==0)cat=203;else if(npr_true==1)cat=204;else cat=205;}
                            else cat=206;
                        }else cat=210;
                    }else cat=998;
                }
            }return cat;
        },{"nu_pdg","ccnc","mc_n_charged_pions_true","mc_n_protons_true","inclusive_strangeness_multiplicity_type"});
    return df_with_event_category.Alias("event_category","event_category_val");
}

ROOT::RDF::RInterface ProcessNuMuVariables(ROOT::RDF::RInterface df, const std::string&){
    auto df_mu_cons = df.Define("muon_mcs_consistency_vec",
        [](const ROOT::RVec<float>& m,const ROOT::RVec<float>& r){ROOT::RVec<float> c(m.size());for(size_t i=0;i<m.size();++i)
        {if(r[i]!=0.f&&std::isfinite(r[i])&&std::isfinite(m[i]))c[i]=(m[i]-r[i])/r[i];else c[i]=-999.f;}return c;},{"trk_mcs_muon_mom_v","trk_range_muon_mom_v"});
    auto df_mu_mask = df_mu_cons.Define("muon_candidate_selection_mask_vec",
        [](const ROOT::RVec<float>& ts,const ROOT::RVec<float>& pid,const ROOT::RVec<float>& sx,const ROOT::RVec<float>& ex,
           const ROOT::RVec<float>& sy,const ROOT::RVec<float>& ey,const ROOT::RVec<float>& sz,const ROOT::RVec<float>& ez,
           const ROOT::RVec<float>& l,const ROOT::RVec<float>& dist,const ROOT::RVec<int>& gen,const ROOT::RVec<float>& mcsc){
            ROOT::RVec<bool> mask(ts.size());for(size_t i=0;i<ts.size();++i){
            bool fid=(sx[i]>5.f&&sx[i]<251.f&&sy[i]>-110.f&&sy[i]<110.f&&sz[i]>20.f&&sz[i]<986.f && ex[i]>5.f&&ex[i]<251.f&&ey[i]>-110.f&&ey[i]<110.f&&ez[i]>20.f&&ez[i]<986.f);
            bool mcok=(mcsc[i]>-0.5f&&mcsc[i]<0.5f);
            mask[i]=(ts[i]>0.8f)&&(pid[i]>0.2f)&&fid&&(l[i]>10.f)&&(dist[i]<4.f)&&(gen[i]==2)&&mcok;}return mask;},
            {"trk_score_v","trk_llr_pid_score_v","trk_sce_start_x_v","trk_sce_end_x_v","trk_sce_start_y_v","trk_sce_end_y_v",
            "trk_sce_start_z_v","trk_sce_end_z_v","trk_len_v","trk_distance_v","pfp_generation_v","muon_mcs_consistency_vec"});
    auto df_sel_idx = df_mu_mask.Define("selected_muon_idx",
        [](const ROOT::RVec<float>&l,const ROOT::RVec<bool>&m){return GetIndexFromVectorSort(l,m,0,false);},{"trk_len_v","muon_candidate_selection_mask_vec"});
    auto df_mu_props = df_sel_idx
    .Define("selected_muon_length",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-1.f);},{"trk_len_v","selected_muon_idx"})
    .Define("selected_muon_momentum_range",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-1.f);},{"trk_range_muon_mom_v","selected_muon_idx"})
    .Define("selected_muon_momentum_mcs",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-1.f);},{"trk_mcs_muon_mom_v","selected_muon_idx"})
    .Define("selected_muon_mcs_consistency",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-999.f);},{"muon_mcs_consistency_vec","selected_muon_idx"})
    .Define("selected_muon_phi",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-999.f);},{"trk_phi_v","selected_muon_idx"})
    .Define("selected_muon_cos_theta",[](const ROOT::RVec<float>&v,int i){float t=GetElementFromVector(v,i,-999.f);return(t>-998.f&&std::isfinite(t))?std::cos(t):-999.f;},{"trk_theta_v","selected_muon_idx"})
    .Define("selected_muon_energy",[](float mom){const float M=0.105658f;return(mom>=0.f&&std::isfinite(mom))?std::sqrt(mom*mom+M*M):-1.f;},{"selected_muon_momentum_range"})
    .Define("selected_muon_start_x",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-999.f);},{"trk_sce_start_x_v","selected_muon_idx"})
    .Define("selected_muon_start_y",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-999.f);},{"trk_sce_start_y_v","selected_muon_idx"})
    .Define("selected_muon_start_z",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-999.f);},{"trk_sce_start_z_v","selected_muon_idx"})
    .Define("selected_muon_end_x",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-999.f);},{"trk_sce_end_x_v","selected_muon_idx"})
    .Define("selected_muon_end_y",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-999.f);},{"trk_sce_end_y_v","selected_muon_idx"})
    .Define("selected_muon_end_z",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-999.f);},{"trk_sce_end_z_v","selected_muon_idx"})
    .Define("selected_muon_trk_score",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-1.f);},{"trk_score_v","selected_muon_idx"})
    .Define("selected_muon_llr_pid_score",[](const ROOT::RVec<float>&v,int i){return GetElementFromVector(v,i,-999.f);},{"trk_llr_pid_score_v","selected_muon_idx"});
    return df_mu_props
    .Define("n_muon_candidates",[](const ROOT::RVec<bool>&m){return ROOT::VecOps::Sum(m);},{"muon_candidate_selection_mask_vec"})
    .Define("n_reco_tracks_passing_basic_cuts",[](const ROOT::RVec<float>&s){int c=0;for(float v:s)if(v>0.1f)c++;return c;},{"trk_score_v"});
}

ROOT::RDF::RInterface AddMCWeightVariables(ROOT::RDF::RInterface df) {
    bool has_nom_w = false; for(const auto&c:df.GetColumnNames()) if(c=="weightSplineTimesTune") has_nom_w=true;
    if(has_nom_w) return df.Define("base_event_weight",[](float w){return(w>0.f&&w<100.f&&std::isfinite(w))?w:1.f;},{"weightSplineTimesTune"});
    return df.Define("base_event_weight",[](){return 1.0f;});
}

ROOT::RDF::RInterface RemoveDuplicates(ROOT::RDF::RInterface df) { return df; }

ROOT::RDF::RInterface LoadSample(
    const ConfigurationManager& config_mgr, const VariableManager& var_mgr,
    const std::string& run_id_key, const std::string& dataset_name_in_code,
    const std::string& variation_name = "cv",
    bool load_reco_event_level = true, bool load_reco_track_level = false, bool load_mc_truth_level = true,
    bool load_mc_weights_and_systematics = true, bool apply_specific_is_event_weights = false,
    const std::vector<std::string>& /* specific_keep_columns - not used if empty */ = {}, bool is_detvar_sample = false
) {
    const auto& paths = config_mgr.GetFilePaths(); std::string prefix, rel_path; NTupleFileType type_enum;
    std::string cat_key = dataset_name_in_code;

    if (is_detvar_sample) {
        const auto& dr_cfg=config_mgr.GetDetVarRunConfig(run_id_key); rel_path=dr_cfg.relative_path_detvar_base;
        const auto& v_map_it=dr_cfg.variations.find(variation_name); if(v_map_it==dr_cfg.variations.end())throw std::runtime_error("V:"+variation_name);
        const auto& d_map_it=v_map_it->second.find(dataset_name_in_code); if(d_map_it==v_map_it->second.end())throw std::runtime_error("D:"+dataset_name_in_code);
        const auto& dv_info=d_map_it->second; const auto& cv_cfg=config_mgr.GetRunConfig(run_id_key);
        const auto& cv_ds_it=cv_cfg.datasets.find(dv_info.name_in_code_cv_equivalent); if(cv_ds_it==cv_cfg.datasets.end())throw std::runtime_error("CV:"+dv_info.name_in_code_cv_equivalent);
        prefix=cv_ds_it->second.file_prefix+dv_info.file_prefix_variant_suffix; type_enum=dv_info.base_sample_type;
        cat_key = IsIntrinsicStrangenessSignalMC(type_enum)?"intrinsic_strangeness_signal":"inclusive_generic_mc";
    } else {
        const auto& r_cfg=config_mgr.GetRunConfig(run_id_key); rel_path=r_cfg.relative_path;
        const auto& ds_it=r_cfg.datasets.find(dataset_name_in_code); if(ds_it==r_cfg.datasets.end())throw std::runtime_error("DS:"+dataset_name_in_code);
        const auto& ds_info=ds_it->second; prefix=ds_info.file_prefix; type_enum=ds_info.type;
        if(IsBeamOnData(type_enum))cat_key="data";else if(IsBeamOffBackground(type_enum))cat_key="ext_background";
        else if(type_enum==NTupleFileType::kMCDirtBNB||type_enum==NTupleFileType::kMCDirtNuMI)cat_key="dirt_mc";
        else if(IsIntrinsicStrangenessSignalMC(type_enum))cat_key="intrinsic_strangeness_signal";
        else if(IsMonteCarlo(type_enum))cat_key="inclusive_generic_mc";else cat_key="unknown_sample_type";
    }
    std::string file_path=paths.ntuple_path_base+"/"+rel_path+"/"+prefix+".root";
    if(!std::ifstream(file_path.c_str()).good())throw std::runtime_error("No file:"+file_path);
    
    ROOT::RDataFrame df_from_file("nuselection/NeutrinoSelectionFilter",file_path);
    ROOT::RDF::RInterface current_df = df_from_file;

    // VariableManager determines which columns are needed by the subsequent processing functions
    // RDataFrame will activate them on demand.
    var_mgr.GetVariables(load_reco_event_level, load_reco_track_level, load_mc_truth_level,
                         load_mc_weights_and_systematics, apply_specific_is_event_weights, type_enum);


    if(load_reco_track_level) current_df=ProcessNuMuVariables(current_df,file_path);
    if(IsMonteCarlo(type_enum)&&!IsDetectorVariation(type_enum)){
        current_df=AddMCWeightVariables(current_df);
        if(apply_specific_is_event_weights&&IsIntrinsicStrangenessSignalMC(type_enum)){
            const std::string is_w_branch="is_channel_correction_weight"; // From VariableManager
            bool has_is_w=false; for(const auto&c:current_df.GetColumnNames())if(c==is_w_branch){has_is_w=true;break;}
            if(has_is_w)current_df=current_df.Redefine("base_event_weight","base_event_weight*"+is_w_branch);
        }
    }
    
    if(IsMonteCarlo(type_enum)||IsBeamOnData(type_enum)||IsBeamOffBackground(type_enum)){
        current_df=AddEventCategories(current_df,cat_key);
    }
    if(IsBeamOffBackground(type_enum)||IsBeamOnData(type_enum))current_df=RemoveDuplicates(current_df);
    return current_df;
}

struct LoadedRunData{std::map<std::string,ROOT::RDF::RInterface>dfs;std::map<std::string,double>scalers;double pot=0.;int trig=0;};

LoadedRunData InternalLoadRun(
    const ConfigurationManager& cfm, const VariableManager& vm, const std::string& rid, const std::string& dkey, const std::string& ekey,
    const std::string& mckey, const std::vector<std::string>& tfkeys, bool blind=true, bool load_is=false, const std::string& iskey="is_bnb",
    bool rl=true, bool rtl=true, bool mctl=true, bool mcws=true, bool apply_isw=true, const std::vector<std::string>& keep_cols={}
){
    LoadedRunData res; const auto& rc=cfm.GetRunConfig(rid);
    const auto& d_info_it=rc.datasets.find(dkey); if(d_info_it==rc.datasets.end())throw std::runtime_error("NoData:"+dkey+"inR:"+rid);
    res.pot=d_info_it->second.pot; res.trig=d_info_it->second.triggers;
    if(!blind) {
        res.dfs[dkey]=LoadSample(cfm,vm,rid,dkey,"cv",rl,rtl,false,false,false,keep_cols,false);
        res.dfs[dkey]=res.dfs[dkey].Define("final_event_weight",[](){return 1.0;});
    }
    res.scalers[dkey]=1.0;
    const auto& e_info_it=rc.datasets.find(ekey); if(e_info_it==rc.datasets.end())throw std::runtime_error("NoEXT:"+ekey+"inR:"+rid);
    bool ext_is_mc=IsMonteCarlo(e_info_it->second.type);
    res.dfs[ekey]=LoadSample(cfm,vm,rid,ekey,"cv",rl,rtl,ext_is_mc?mctl:false,ext_is_mc?mcws:false,false,keep_cols,false);
    double es=1.; if(e_info_it->second.triggers>0&&res.trig>0)es=(double)res.trig/e_info_it->second.triggers;
    else if(ext_is_mc&&e_info_it->second.pot>0&&res.pot>0)es=res.pot/e_info_it->second.pot;
    if(ext_is_mc)res.dfs[ekey]=res.dfs[ekey].Define("final_event_weight",[es](float bew){return bew*es;},{"base_event_weight"});
    else res.dfs[ekey]=res.dfs[ekey].Define("final_event_weight",[es](){return es;});
    res.scalers[ekey]=es;

    std::map<std::string,ROOT::RDF::RInterface> tmp_mc_dfs; std::string main_mc_excl_filt="1==1";
    std::vector<std::string> all_mc_keys={mckey}; all_mc_keys.insert(all_mc_keys.end(),tfkeys.begin(),tfkeys.end());
    if(load_is)all_mc_keys.push_back(iskey);
    std::sort(all_mc_keys.begin(),all_mc_keys.end());all_mc_keys.erase(std::unique(all_mc_keys.begin(),all_mc_keys.end()),all_mc_keys.end());

    for(const auto&mk:all_mc_keys){
        if(mk.empty())continue; const auto&mc_ds_it=rc.datasets.find(mk); if(mc_ds_it==rc.datasets.end())throw std::runtime_error("NoMC:"+mk+"inR:"+rid);
        const auto&mc_ds_info=mc_ds_it->second; bool is_this_is=(mk==iskey&&load_is);
        tmp_mc_dfs[mk]=LoadSample(cfm,vm,rid,mk,"cv",rl,rtl,mctl,mcws,apply_isw&&is_this_is,keep_cols,false);
        double mc_s=(mc_ds_info.pot>0&&res.pot>0)?res.pot/mc_ds_info.pot:1.; res.scalers[mk]=mc_s;
        tmp_mc_dfs[mk]=tmp_mc_dfs[mk].Define("final_event_weight",[mc_s](float b){return b*mc_s;},{"base_event_weight"});
        if(mk!=mckey&&!mc_ds_info.filter_query.empty()){if(main_mc_excl_filt=="1==1")main_mc_excl_filt="!("+mc_ds_info.filter_query+")";else main_mc_excl_filt+="&&!("+mc_ds_info.filter_query+")";}
    }
    for(const auto&k:tfkeys)if(tmp_mc_dfs.count(k))res.dfs[k]=tmp_mc_dfs[k];
    if(load_is&&tmp_mc_dfs.count(iskey))res.dfs[iskey]=tmp_mc_dfs[iskey];
    if(tmp_mc_dfs.count(mckey)){if(main_mc_excl_filt!="1==1")res.dfs[mckey]=tmp_mc_dfs[mckey].Filter(main_mc_excl_filt,"MainMCExcl");else res.dfs[mckey]=tmp_mc_dfs[mckey];}
    return res;
}

LoadedRunData LoadRuns(
    const ConfigurationManager& cfm, const VariableManager& vm, const std::vector<std::string>& rids,
    const std::string& d_pat="data_%B%", const std::string& e_pat="ext_%B%", const std::string& mc_pat="mc_%B%",
    const std::vector<std::string>& tf_pats={}, bool blind=true, bool load_is=false, const std::string& is_pat="is_%B%",
    bool rl=true, bool rtl=true, bool mctl=true, bool mcws=true, bool apply_isw=true
){
    LoadedRunData comb_res; double tot_pot=0.; std::map<std::string,std::vector<ROOT::RDF::RInterface>> df_lists;
    std::map<std::string,std::vector<double>> mc_pots_map;
    for(const auto&rid:rids){
        const auto&rc=cfm.GetRunConfig(rid); std::string bt=(rc.beam_type=="BNB")?"bnb":"numi";
        auto rep=[&](std::string p){size_t pos=p.find("%B%");if(pos!=std::string::npos)p.replace(pos,3,bt);return p;};
        std::vector<std::string> cur_tf_keys; for(const auto&p:tf_pats)cur_tf_keys.push_back(rep(p));
        LoadedRunData srd=InternalLoadRun(cfm,vm,rid,rep(d_pat),rep(e_pat),rep(mc_pat),cur_tf_keys,blind,load_is,rep(is_pat),rl,rtl,mctl,mcws,apply_isw);
        tot_pot+=srd.pot; for(auto const&[k,df]:srd.dfs)if(df.GetNRuns()>0)df_lists[k].push_back(df);
        for(auto const&[k,s]:srd.scalers){
            if(k!=rep(d_pat)&&k!=rep(e_pat)){const auto&dsit=rc.datasets.find(k);if(dsit!=rc.datasets.end())mc_pots_map[k].push_back(dsit->second.pot);}
            else comb_res.scalers[k]=s;
        }
    }
    for(auto&p:df_lists){if(!p.second.empty())comb_res.dfs[p.first]=(p.second.size()>1)?ROOT::RDF::Experimental::Concatenate(p.second):p.second[0];}
    comb_res.pot=tot_pot;
    for(auto const&[k,pv]:mc_pots_map){double totmc_p=std::accumulate(pv.begin(),pv.end(),0.);comb_res.scalers[k]=(totmc_p>0&&tot_pot>0)?tot_pot/totmc_p:1.;}
    return comb_res;
}

} // namespace DataLoader
} // namespace AnalysisFramework


