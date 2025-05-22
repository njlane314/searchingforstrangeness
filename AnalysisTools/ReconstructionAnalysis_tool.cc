#ifndef ANALYSIS_PANDORATRACKANALYSIS_CXX
#define ANALYSIS_PANDORATRACKANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "../CommonDefs/Types.h" 
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/LLR_PID.h"
#include "../CommonDefs/LLRPID_proton_muon_lookup.h"
#include "../CommonDefs/LLRPID_correction_lookup.h"
#include "../CommonDefs/CalibrationFuncs.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "TTree.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>
#include <optional> 

namespace analysis 
{
    class ReconstructionAnalysis : public AnalysisToolBase {
    public:
        explicit ReconstructionAnalysis(fhicl::ParameterSet const& p);
        virtual ~ReconstructionAnalysis() = default;

        void configure(const fhicl::ParameterSet& p) override;
        void analyseEvent(art::Event const& e, bool fData) override; 
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override; 
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;

    private:
        art::InputTag _fPFPproducer;
        art::InputTag _fTRKproducer;
        art::InputTag _fCALOproducer;
        art::InputTag _fHproducer;
        art::InputTag _fBacktrackTag;
        art::InputTag _fMCRproducer;

        bool _fRecalibrateHits;
        float _fEnergyThresholdForMCHits;
        std::vector<float> _fADCtoE;

        trkf::TrackMomentumCalculator _trkmom;
        trkf::TrajectoryMCSFitter _mcsfitter;
        TParticlePDG* _proton_pdg;
        TParticlePDG* _muon_pdg;
        common::LLRPID _llr_pid_calculator;
        common::ProtonMuonLookUpParameters _protonmuon_parameters;
        common::CorrectionLookUpParameters _correction_parameters;

        float _reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z;
        float _topological_score;

        std::vector<size_t> _trk_pfp_id_v;
        std::vector<float> _trk_score_v;
        std::vector<float> _trk_len_v;
        std::vector<float> _trk_distance_v;
        std::vector<float> _trk_theta_v;
        std::vector<float> _trk_phi_v;
        std::vector<float> _trk_llr_pid_score_v;
        std::vector<float> _trk_mcs_muon_mom_v;
        std::vector<float> _trk_range_muon_mom_v;
        std::vector<int> _trk_generation_v;
        std::vector<float> _trk_llr_pid_sum_v; 

        void fillDefault(); 
    };

    ReconstructionAnalysis::ReconstructionAnalysis(const fhicl::ParameterSet& pset) :
        _mcsfitter(fhicl::Table<trkf::TrajectoryMCSFitter::Config>(pset.get<fhicl::ParameterSet>("mcsfitmu"))) 
    {
        this->configure(pset);

        _proton_pdg = TDatabasePDG::Instance()->GetParticle(2212);
        _muon_pdg = TDatabasePDG::Instance()->GetParticle(13);

        _llr_pid_calculator.set_dedx_binning(0, _protonmuon_parameters.dedx_edges_pl_0);
        _llr_pid_calculator.set_par_binning(0, _protonmuon_parameters.parameters_edges_pl_0);
        _llr_pid_calculator.set_lookup_tables(0, _protonmuon_parameters.dedx_pdf_pl_0);
        _llr_pid_calculator.set_dedx_binning(1, _protonmuon_parameters.dedx_edges_pl_1);
        _llr_pid_calculator.set_par_binning(1, _protonmuon_parameters.parameters_edges_pl_1);
        _llr_pid_calculator.set_lookup_tables(1, _protonmuon_parameters.dedx_pdf_pl_1);
        _llr_pid_calculator.set_dedx_binning(2, _protonmuon_parameters.dedx_edges_pl_2);
        _llr_pid_calculator.set_par_binning(2, _protonmuon_parameters.parameters_edges_pl_2);
        _llr_pid_calculator.set_lookup_tables(2, _protonmuon_parameters.dedx_pdf_pl_2);
        
        if (_fRecalibrateHits) {
            _llr_pid_calculator.set_corr_par_binning(0, _correction_parameters.parameter_correction_edges_pl_0);
            _llr_pid_calculator.set_correction_tables(0, _correction_parameters.correction_table_pl_0);
            _llr_pid_calculator.set_corr_par_binning(1, _correction_parameters.parameter_correction_edges_pl_1);
            _llr_pid_calculator.set_correction_tables(1, _correction_parameters.correction_table_pl_1);
            _llr_pid_calculator.set_corr_par_binning(2, _correction_parameters.parameter_correction_edges_pl_2);
            _llr_pid_calculator.set_correction_tables(2, _correction_parameters.correction_table_pl_2);
        }
    }

    void ReconstructionAnalysis::configure(const fhicl::ParameterSet& p) {
        _fPFPproducer = p.get<art::InputTag>("PFPproducer", "pandora");
        _fTRKproducer = p.get<art::InputTag>("TRKproducer", "pandora");
        _fCALOproducer = p.get<art::InputTag>("CALOproducer", "pandoracalo");
        
        _fRecalibrateHits = p.get<bool>("RecalibrateHits", false);
        if (_fRecalibrateHits) {
            _fHproducer = p.get<art::InputTag>("Hproducer", "gaushit");
            _fBacktrackTag = p.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
            _fMCRproducer = p.get<art::InputTag>("MCRproducer", "mcreco");
            _fEnergyThresholdForMCHits = p.get<float>("EnergyThresholdForMCHits", 0.1);
        }
        _fADCtoE = p.get<std::vector<float>>("ADCtoE");
    }

    void ReconstructionAnalysis::analyseEvent(art::Event const& e, bool fData) {}

    void ReconstructionAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) {
        _reco_nu_vtx_x = std::numeric_limits<float>::lowest();
        _reco_nu_vtx_y = std::numeric_limits<float>::lowest();
        _reco_nu_vtx_z = std::numeric_limits<float>::lowest();
        _topological_score = std::numeric_limits<float>::lowest();
        
        TVector3 nuvtx_internal(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        size_t primary_neutrino_pfp_id = std::numeric_limits<size_t>::max();

        for (const auto& pfp_elem : slice_pfp_v) {
            if (pfp_elem->IsPrimary() && (std::abs(pfp_elem->PdgCode()) == 12 || std::abs(pfp_elem->PdgCode()) == 14 || std::abs(pfp_elem->PdgCode()) == 16)) {
                primary_neutrino_pfp_id = pfp_elem->Self();

                auto vtx_assns_view = pfp_elem.get<recob::Vertex>(); 
                if (vtx_assns_view.size() > 0) {
                    art::Ptr<recob::Vertex> vtx_ptr = vtx_assns_view.at(0); 
                    if (vtx_ptr) { 
                        double xyz[3] = {0.0, 0.0, 0.0};
                        vtx_ptr->XYZ(xyz);
                        
                        float sce_corrected_xyz[3];
                        common::ApplySCECorrectionXYZ(xyz[0], xyz[1], xyz[2], sce_corrected_xyz);
                        _reco_nu_vtx_x = sce_corrected_xyz[0];
                        _reco_nu_vtx_y = sce_corrected_xyz[1];
                        _reco_nu_vtx_z = sce_corrected_xyz[2];
                        nuvtx_internal.SetXYZ(_reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z);
                    }
                }

                auto metadata_assns_view = pfp_elem.get<larpandoraobj::PFParticleMetadata>();
                if (metadata_assns_view.size() > 0) {
                    art::Ptr<larpandoraobj::PFParticleMetadata> metadata_ptr = metadata_assns_view.at(0); 
                    if (metadata_ptr) { 
                        auto pfParticlePropertiesMap = metadata_ptr->GetPropertiesMap();
                        if (pfParticlePropertiesMap.count("NuScore")) {
                            _topological_score = pfParticlePropertiesMap.at("NuScore");
                        }
                    }
                }
                break;
            }
        }
        
        common::ProxyCaloColl_t const* calo_proxy_ptr = nullptr;
        std::optional<common::ProxyCaloColl_t> calo_proxy_holder;
        bool calo_proxy_ok = false;
        try {
            calo_proxy_holder.emplace(proxy::getCollection<std::vector<recob::Track>>(e, _fTRKproducer, proxy::withAssociated<anab::Calorimetry>(_fCALOproducer)));
            calo_proxy_ptr = &calo_proxy_holder.value();
            calo_proxy_ok = true;
        } catch (const cet::exception& ex) {
            mf::LogWarning("ReconstructionAnalysis") << "Could not get Track/Calo collections: " << ex.what();
        }
        
        std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
        if (!fData && _fRecalibrateHits) {
            art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(_fHproducer);
            assocMCPart = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(inputHits, e, _fBacktrackTag);
        }

        for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); ++i_pfp) {
            const auto& pfp = slice_pfp_v[i_pfp]; 
            if (pfp->IsPrimary()) continue;

            auto trk_assns_view = pfp.get<recob::Track>();
            if (trk_assns_view.size() == 1) {
                art::Ptr<recob::Track> trk = trk_assns_view.at(0); 
                if (trk && calo_proxy_ok) { 
                    _trk_pfp_id_v.push_back(pfp->Self());

                    _trk_len_v.push_back(common::GetSCECorrTrackLength(trk));
                    _trk_theta_v.push_back(trk->Theta());
                    _trk_phi_v.push_back(trk->Phi());

                    TVector3 trk_start_v3_sce(0,0,0);
                    float start_sce[3]; 
                    common::ApplySCECorrectionXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z(), start_sce);
                    trk_start_v3_sce.SetXYZ(start_sce[0], start_sce[1], start_sce[2]);

                    if (!std::isnan(nuvtx_internal.X())) {
                        _trk_distance_v.push_back((trk_start_v3_sce - nuvtx_internal).Mag());
                    } else {
                        _trk_distance_v.push_back(std::numeric_limits<float>::lowest());
                    }

                    float mcs_mom = -1.f, range_mom_mu = -1.f;
                    if (_muon_pdg) {
                        mcs_mom = _mcsfitter.fitMcs(trk->Trajectory(), _muon_pdg->PdgCode()).bestMomentum();
                        range_mom_mu = _trkmom.GetTrackMomentum(common::GetSCECorrTrackLength(trk), _muon_pdg->PdgCode());
                    }
                    _trk_mcs_muon_mom_v.push_back(mcs_mom);
                    _trk_range_muon_mom_v.push_back(range_mom_mu);

                    float llr_pid_sum = 0.f;
                    auto calo_data_view = (*calo_proxy_ptr)[trk.key()].get<anab::Calorimetry>();
                    for (const art::Ptr<anab::Calorimetry>& calo_ptr : calo_data_view) { 
                        if (!calo_ptr) continue; 
                        auto const& plane_id = calo_ptr->PlaneID().Plane;
                        if (plane_id < 0 || static_cast<size_t>(plane_id) >= _fADCtoE.size()) continue;

                        std::vector<float> dqdx_corr = calo_ptr->dQdx();
                        if (!fData && _fRecalibrateHits && assocMCPart) {
                            dqdx_corr = _llr_pid_calculator.correct_many_hits_one_plane(calo_ptr, *trk, assocMCPart, _fRecalibrateHits, _fEnergyThresholdForMCHits, false);
                        }
                        
                        std::vector<float> dedx_corr;
                        for (size_t i_hit = 0; i_hit < dqdx_corr.size(); ++i_hit) {
                            if (i_hit >= calo_ptr->XYZ().size()) continue;
                            float dedx_val = common::ModBoxCorrection(dqdx_corr[i_hit] * _fADCtoE[plane_id], calo_ptr->XYZ()[i_hit].X(), calo_ptr->XYZ()[i_hit].Y(), calo_ptr->XYZ()[i_hit].Z());
                            dedx_corr.push_back(dedx_val);
                        }

                        std::vector<std::vector<float>> par_values_for_llr;
                        par_values_for_llr.push_back(calo_ptr->ResidualRange());
                        par_values_for_llr.push_back(calo_ptr->TrkPitchVec());
                        
                        if (!dedx_corr.empty() && dedx_corr.size() == calo_ptr->ResidualRange().size() && dedx_corr.size() == calo_ptr->TrkPitchVec().size()){
                            llr_pid_sum += _llr_pid_calculator.LLR_many_hits_one_plane(dedx_corr, par_values_for_llr, plane_id);
                        }
                    }
                    _trk_llr_pid_sum_v.push_back(llr_pid_sum);
                    _trk_llr_pid_score_v.push_back(std::atan(llr_pid_sum / 100.f) * 2.f / M_PI);

                    float current_trk_score = std::numeric_limits<float>::lowest();
                    auto trk_metadata_assns_view = pfp.get<larpandoraobj::PFParticleMetadata>();
                    if (trk_metadata_assns_view.size() > 0) {
                        art::Ptr<larpandoraobj::PFParticleMetadata> trk_metadata_ptr = trk_metadata_assns_view.at(0); 
                        if (trk_metadata_ptr) { 
                            auto trkParticlePropertiesMap = trk_metadata_ptr->GetPropertiesMap();
                            if (trkParticlePropertiesMap.count("TrackScore")) {
                                current_trk_score = trkParticlePropertiesMap.at("TrackScore");
                            } else if (trkParticlePropertiesMap.count("TaggerScore")) {
                                current_trk_score = trkParticlePropertiesMap.at("TaggerScore");
                            }
                        }
                    }
                    _trk_score_v.push_back(current_trk_score);

                    int generation = 0;
                    if (primary_neutrino_pfp_id != std::numeric_limits<size_t>::max()) {
                        if (pfp->Parent() == primary_neutrino_pfp_id) {
                            generation = 2;
                        }
                    }
                    _trk_generation_v.push_back(generation);

                } else { // if !trk or !calo_proxy_ok
                    this->fillDefault();
                }
            } else { // if trk_assns_view.size() != 1
                this->fillDefault();
            }
        }
    }

    void ReconstructionAnalysis::fillDefault() {
        _trk_pfp_id_v.push_back(std::numeric_limits<size_t>::max());
        _trk_score_v.push_back(std::numeric_limits<float>::lowest());
        _trk_len_v.push_back(std::numeric_limits<float>::lowest());
        _trk_distance_v.push_back(std::numeric_limits<float>::lowest());
        _trk_theta_v.push_back(std::numeric_limits<float>::lowest());
        _trk_phi_v.push_back(std::numeric_limits<float>::lowest());
        _trk_llr_pid_sum_v.push_back(std::numeric_limits<float>::lowest());
        _trk_llr_pid_score_v.push_back(std::numeric_limits<float>::lowest());
        _trk_mcs_muon_mom_v.push_back(std::numeric_limits<float>::lowest());
        _trk_range_muon_mom_v.push_back(std::numeric_limits<float>::lowest());
        _trk_generation_v.push_back(std::numeric_limits<int>::lowest());
    }

    void ReconstructionAnalysis::setBranches(TTree* _tree) {
        if (!_tree) return;
        _tree->Branch("reco_nu_vtx_x", &_reco_nu_vtx_x, "reco_nu_vtx_x/F");
        _tree->Branch("reco_nu_vtx_y", &_reco_nu_vtx_y, "reco_nu_vtx_y/F");
        _tree->Branch("reco_nu_vtx_z", &_reco_nu_vtx_z, "reco_nu_vtx_z/F");
        _tree->Branch("topological_score", &_topological_score, "topological_score/F");

        _tree->Branch("trk_pfp_id_v", &_trk_pfp_id_v);
        _tree->Branch("trk_score_v", &_trk_score_v);
        _tree->Branch("trk_len_v", &_trk_len_v);
        _tree->Branch("trk_distance_v", &_trk_distance_v);
        _tree->Branch("trk_theta_v", &_trk_theta_v);
        _tree->Branch("trk_phi_v", &_trk_phi_v);
        _tree->Branch("trk_llr_pid_sum_v", &_trk_llr_pid_sum_v);
        _tree->Branch("trk_llr_pid_score_v", &_trk_llr_pid_score_v);
        _tree->Branch("trk_mcs_muon_mom_v", &_trk_mcs_muon_mom_v);
        _tree->Branch("trk_range_muon_mom_v", &_trk_range_muon_mom_v);
        _tree->Branch("trk_generation_v", &_trk_generation_v);
    }

    void ReconstructionAnalysis::resetTTree(TTree* _tree) {
        _reco_nu_vtx_x = std::numeric_limits<float>::lowest(); 
        _reco_nu_vtx_y = std::numeric_limits<float>::lowest(); 
        _reco_nu_vtx_z = std::numeric_limits<float>::lowest();
        _topological_score = std::numeric_limits<float>::lowest();

        _trk_pfp_id_v.clear();
        _trk_score_v.clear();
        _trk_len_v.clear();
        _trk_distance_v.clear();
        _trk_theta_v.clear();
        _trk_phi_v.clear();
        _trk_llr_pid_sum_v.clear();
        _trk_llr_pid_score_v.clear();
        _trk_mcs_muon_mom_v.clear();
        _trk_range_muon_mom_v.clear();
        _trk_generation_v.clear();
    }

    DEFINE_ART_CLASS_TOOL(ReconstructionAnalysis)
} 
#endif

