#ifndef SELECTION_NUMUSELECTION_CXX
#define SELECTION_NUMUSELECTION_CXX

#include "SelectionToolBase.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/Types.h"
#include "larcore/Geometry/Geometry.h"
#include "../CommonDefs/LLR_PID.h"
#include "../CommonDefs/LLRPID_proton_muon_lookup.h"
#include "../CommonDefs/LLRPID_correction_lookup.h"
#include "../CommonDefs/CalibrationFuncs.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

namespace selection 
{
    class PreSelection : public SelectionToolBase {
    public:
        PreSelection(const fhicl::ParameterSet &pset);
        ~PreSelection();

        void configure(const fhicl::ParameterSet &pset);

        bool selectEvent(art::Event const &e, const std::vector<common::ProxyPfpElem_t> &pfp_pxy_v) override;
        
        void setBranches(TTree *_tree) override;
        
        void resetTTree(TTree *_tree) override;

    private:
        bool sel_vertex_in_fv;
        bool sel_pfp_track_passed;
        bool sel_nu_mu_cc;
        bool sel_has_mu_candidate;
        bool sel_topo_score_passed;
        int mu_candidate_trk_key;
        float slice_nu_score;
        float fTopoScoreCut;

        art::InputTag fTRKproducer;
        art::InputTag fPIDproducer;
        art::InputTag fCALOproducer;
        art::InputTag fHproducer;
        art::InputTag fBacktrackTag;
        float fTrkScoreCut;
        float fVertexDistanceCut;
        float fTrackLengthCut;
        float fLLRPIDScoreCut;
        bool fRecalibrateHits;
        float fEnergyThresholdForMCHits;
        std::vector<float> fADCtoE;
        
        float fFidvolXstart;
        float fFidvolXend;
        float fFidvolYstart;
        float fFidvolYend;
        float fFidvolZstart;
        float fFidvolZend;

        common::LLRPID llr_pid_calculator;
        common::ProtonMuonLookUpParameters protonmuon_parameters;
        common::CorrectionLookUpParameters correction_parameters;

        bool isFiducial(const double x[3]) const;
    };

    PreSelection::PreSelection(const fhicl::ParameterSet &pset) {
        this->configure(pset);
    }

    PreSelection::~PreSelection() {}

    void PreSelection::configure(const fhicl::ParameterSet &pset) {
        fTRKproducer = pset.get<art::InputTag>("TRKproducer", "pandoraTrack");
        fPIDproducer = pset.get<art::InputTag>("PIDproducer", "pandoraTrackcalipid");
        fCALOproducer = pset.get<art::InputTag>("CALOproducer", "pandoraTrackcali");
        fHproducer = pset.get<art::InputTag>("Hproducer", "gaushit");
        fBacktrackTag = pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
        fTrkScoreCut = pset.get<float>("TrkScoreCut", 0.8);
        fVertexDistanceCut = pset.get<float>("VertexDistanceCut", 4.0);
        fTrackLengthCut = pset.get<float>("TrackLengthCut", 10.0);
        fLLRPIDScoreCut = pset.get<float>("LLRPIDScoreCut", 0.2);
        fRecalibrateHits = pset.get<bool>("RecalibrateHits", false);
        fEnergyThresholdForMCHits = pset.get<float>("EnergyThresholdForMCHits", 0.1);
        fADCtoE = pset.get<std::vector<float>>("ADCtoE");
        fTopoScoreCut = pset.get<float>("TopoScoreCut", 0.2);

        fFidvolXstart = pset.get<float>("FidvolXstart", 10.0);
        fFidvolXend = pset.get<float>("FidvolXend", 10.0);
        fFidvolYstart = pset.get<float>("FidvolYstart", 10.0);
        fFidvolYend = pset.get<float>("FidvolYend", 10.0);
        fFidvolZstart = pset.get<float>("FidvolZstart", 10.0);
        fFidvolZend = pset.get<float>("FidvolZend", 50.0);

        llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
        llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
        llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

        llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
        llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
        llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

        llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
        llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
        llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);

        if (fRecalibrateHits) {
            llr_pid_calculator.set_corr_par_binning(0, correction_parameters.parameter_correction_edges_pl_0);
            llr_pid_calculator.set_correction_tables(0, correction_parameters.correction_table_pl_0);

            llr_pid_calculator.set_corr_par_binning(1, correction_parameters.parameter_correction_edges_pl_1);
            llr_pid_calculator.set_correction_tables(1, correction_parameters.correction_table_pl_1);

            llr_pid_calculator.set_corr_par_binning(2, correction_parameters.parameter_correction_edges_pl_2);
            llr_pid_calculator.set_correction_tables(2, correction_parameters.correction_table_pl_2);
        }
    }

    bool PreSelection::selectEvent(art::Event const &e, const std::vector<common::ProxyPfpElem_t> &pfp_pxy_v) {
        common::ProxyPIDColl_t const &pid_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer, proxy::withAssociated<anab::ParticleID>(fPIDproducer));
        common::ProxyCaloColl_t const &calo_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer, proxy::withAssociated<anab::Calorimetry>(fCALOproducer));

        std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
        if (!e.isRealData() && fRecalibrateHits) {
            art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
            assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
        }

        double nu_vtx[3] = {};
        bool has_neutrino = false;
        size_t neutrino_pfp_index = 0;
        for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++) {
            auto pdg = std::abs(pfp_pxy_v[i_pfp]->PdgCode());
            if (pdg == 12 || pdg == 14) {
                auto vtx_v = pfp_pxy_v[i_pfp].get<recob::Vertex>();
                if (vtx_v.size() == 1) {
                    vtx_v[0]->XYZ(nu_vtx);
                    has_neutrino = true;
                    neutrino_pfp_index = i_pfp;
                    break;
                }
            }
        }
        if (!has_neutrino) {
            slice_nu_score = std::numeric_limits<float>::lowest();
            sel_topo_score_passed = false;
            return false;
        }

        auto metadata_pxy_v = pfp_pxy_v[neutrino_pfp_index].get<larpandoraobj::PFParticleMetadata>();
        if (metadata_pxy_v.size() == 0) {
            slice_nu_score = std::numeric_limits<float>::lowest();
            sel_topo_score_passed = false;
            return false;
        }
        const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(metadata_pxy_v.at(0));
        auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
        if (pfParticlePropertiesMap.find("NuScore") == pfParticlePropertiesMap.end()) {
            slice_nu_score = std::numeric_limits<float>::lowest();
            sel_topo_score_passed = false;
            return false;
        }
        slice_nu_score = pfParticlePropertiesMap["NuScore"];
        if (slice_nu_score <= fTopoScoreCut) {
            sel_topo_score_passed = false;
            return false;
        } else {
            sel_topo_score_passed = true;
        }

        sel_vertex_in_fv = isFiducial(nu_vtx);
        if (!sel_vertex_in_fv) return false;

        bool has_pfp_track = false;
        std::vector<std::pair<size_t, float>> mu_candidates;
        for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++) {
            const auto &pfp_pxy = pfp_pxy_v[i_pfp];
            float trk_score = common::GetTrackShowerScore(pfp_pxy);
            if (trk_score > fTrkScoreCut) {
                auto trk_v = pfp_pxy.get<recob::Track>();
                if (trk_v.size() > 0) {
                    const auto &trk = trk_v[0];
                    double trk_start[3] = {trk->Start().X(), trk->Start().Y(), trk->Start().Z()};
                    double distance = std::sqrt(std::pow(nu_vtx[0] - trk_start[0], 2) +
                                                std::pow(nu_vtx[1] - trk_start[1], 2) +
                                                std::pow(nu_vtx[2] - trk_start[2], 2));
                    if (distance < fVertexDistanceCut) has_pfp_track = true;

                    float trk_len = trk->Length();
                    if (trk_len > fTrackLengthCut) {
                        auto calopxy = calo_proxy[trk.key()];
                        auto calo_v = calopxy.get<anab::Calorimetry>();
                        float llr_pid_total = 0;
                        for (auto const &calo : calo_v) {
                            auto const &plane = calo->PlaneID().Plane;
                            auto const &dqdx_values = calo->dQdx();
                            auto const &rr = calo->ResidualRange();
                            auto const &pitch = calo->TrkPitchVec();
                            auto const &xyz_v = calo->XYZ();
                            std::vector<std::vector<float>> par_values;
                            par_values.push_back(rr);
                            par_values.push_back(pitch);

                            std::vector<float> dqdx_values_corrected, dedx_values_corrected;
                            if (!e.isRealData() && fRecalibrateHits && assocMCPart) {
                                dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(calo, trk.value(), assocMCPart, fRecalibrateHits, fEnergyThresholdForMCHits, false);
                            } else {
                                dqdx_values_corrected = dqdx_values;
                            }

                            for (size_t i = 0; i < dqdx_values_corrected.size(); i++) {
                                float aux_dedx = common::ModBoxCorrection(dqdx_values_corrected[i] * fADCtoE[plane], xyz_v[i].X(), xyz_v[i].Y(), xyz_v[i].Z());
                                dedx_values_corrected.push_back(aux_dedx);
                            }

                            float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);
                            llr_pid_total += llr_pid;
                        }
                        float llr_pid_score = atan(llr_pid_total / 100.) * 2 / 3.14159266;
                        if (llr_pid_score > fLLRPIDScoreCut) {
                            mu_candidates.push_back(std::make_pair(i_pfp, llr_pid_score));
                        }
                    }
                }
            }
        }

        sel_pfp_track_passed = has_pfp_track;
        sel_has_mu_candidate = !mu_candidates.empty();

        if (sel_has_mu_candidate) {
            auto max_it = std::max_element(mu_candidates.begin(), mu_candidates.end(),
                                           [](const auto &a, const auto &b) { return a.second < b.second; });
            size_t muon_index = max_it->first;
            const auto &muon_pfp = pfp_pxy_v[muon_index];
            auto muon_trk_v = muon_pfp.get<recob::Track>();
            if (muon_trk_v.size() > 0) {
                const auto &muon_trk = muon_trk_v[0];
                mu_candidate_trk_key = muon_trk.key();
            }
        } else {
            mu_candidate_trk_key = -1;
        }

        sel_nu_mu_cc = sel_vertex_in_fv && sel_pfp_track_passed && sel_has_mu_candidate;
        return sel_nu_mu_cc;
    }

    void PreSelection::setBranches(TTree *_tree) {
        _tree->Branch("sel_vertex_in_fv", &sel_vertex_in_fv, "sel_vertex_in_fv/O");
        _tree->Branch("sel_pfp_track_passed", &sel_pfp_track_passed, "sel_pfp_track_passed/O");
        _tree->Branch("sel_has_mu_candidate", &sel_has_mu_candidate, "sel_has_mu_candidate/O");
        _tree->Branch("sel_nu_mu_cc", &sel_nu_mu_cc, "sel_nu_mu_cc/O");
        _tree->Branch("mu_candidate_trk_key", &mu_candidate_trk_key, "mu_candidate_trk_key/I");
        _tree->Branch("slice_nu_score", &slice_nu_score, "slice_nu_score/F");
        _tree->Branch("sel_topo_score_passed", &sel_topo_score_passed, "sel_topo_score_passed/O");
    }

    void PreSelection::resetTTree(TTree *_tree) {
        sel_vertex_in_fv = false;
        sel_pfp_track_passed = false;
        sel_has_mu_candidate = false;
        sel_nu_mu_cc = false;
        mu_candidate_trk_key = -1;
        slice_nu_score = std::numeric_limits<float>::lowest();
        sel_topo_score_passed = false;
    }

    bool PreSelection::isFiducial(const double x[3]) const {
        auto const& tpc = art::ServiceHandle<geo::Geometry>{}->TPC();
        std::vector<double> bnd = {0., 2. * tpc.HalfWidth(), -tpc.HalfHeight(), tpc.HalfHeight(), 0., tpc.Length()};
        bool is_x = x[0] > (bnd[0] + fFidvolXstart) && x[0] < (bnd[1] - fFidvolXend);
        bool is_y = x[1] > (bnd[2] + fFidvolYstart) && x[1] < (bnd[3] - fFidvolYend);
        bool is_z = x[2] > (bnd[4] + fFidvolZstart) && x[2] < (bnd[5] - fFidvolZend);
        return is_x && is_y && is_z;
    }

    DEFINE_ART_CLASS_TOOL(PreSelection)
} 

#endif