#ifndef ANALYSIS_SHOWERANALYSIS_CXX
#define ANALYSIS_SHOWERANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/ProximityClustering.h"
#include "../CommonDefs/TrackFitterFunctions.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/LLR_PID.h"
#include "../CommonDefs/LLRPID_correction_lookup.h"
#include "../CommonDefs/LLRPID_electron_photon_lookup.h"
#include "../CommonDefs/Types.h"

namespace analysis
{

class ShowerAnalysis : public AnalysisToolBase
{

public:

    ShowerAnalysis(const fhicl::ParameterSet &pset);
    ~ShowerAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);

    void analyseEvent(art::Event const &e, bool fData) override;

    void analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

    void SaveTruth(art::Event const &e);

    void fillDefault();

    void setBranches(TTree *_tree) override;

    void resetTTree(TTree *_tree) override;

private:
    int _run, _sub, _evt;

    art::InputTag fTRKproducer;
    art::InputTag fCALproducer;

    art::InputTag fBacktrackTag;
    art::InputTag fHproducer;
    float fEnergyThresholdForMCHits;

    float fdEdxcmSkip, fdEdxcmLen;
    bool fLocaldEdx;

    std::vector<float> fADCtoE;
    bool fRecalibrateHits;

    common::LLRPID llr_pid_calculator;
    common::ElectronPhotonLookUpParameters electronphoton_parameters;
    common::CorrectionLookUpParameters correction_parameters;

    std::vector<float> _shr_energy_u_v;
    std::vector<float> _shr_energy_v_v;
    std::vector<float> _shr_energy_y_v;

    std::vector<float> _shr_dedx_u_v;
    std::vector<float> _shr_dedx_v_v;
    std::vector<float> _shr_dedx_y_v;

    std::vector<size_t> _shr_pfp_id_v;

    std::vector<float> _shr_start_x_v;
    std::vector<float> _shr_start_y_v;
    std::vector<float> _shr_start_z_v;

    std::vector<float> _shr_start_U_v;
    std::vector<float> _shr_start_V_v;
    std::vector<float> _shr_dist_v;

    std::vector<float> _shr_px_v;
    std::vector<float> _shr_py_v;
    std::vector<float> _shr_pz_v;

    std::vector<float> _shr_theta_v;
    std::vector<float> _shr_phi_v;

    std::vector<float> _shr_pitch_u_v;
    std::vector<float> _shr_pitch_v_v;
    std::vector<float> _shr_pitch_y_v;

    std::vector<float> _shr_openangle_v;

    std::vector<int> _shr_tkfit_nhits_v;
    std::vector<float> _shr_tkfit_start_x_v;
    std::vector<float> _shr_tkfit_start_y_v;
    std::vector<float> _shr_tkfit_start_z_v;

    std::vector<float> _shr_tkfit_start_U_v;
    std::vector<float> _shr_tkfit_start_V_v;

    std::vector<float> _shr_tkfit_theta_v;
    std::vector<float> _shr_tkfit_phi_v;

    std::vector<float> _shr_tkfit_pitch_u_v;
    std::vector<float> _shr_tkfit_pitch_v_v;
    std::vector<float> _shr_tkfit_pitch_y_v;

    std::vector<float> _shr_tkfit_dedx_u_v;
    std::vector<float> _shr_tkfit_dedx_v_v;
    std::vector<float> _shr_tkfit_dedx_y_v;

    std::vector<float> _shr_tkfit_gap10_dedx_u_v;
    std::vector<float> _shr_tkfit_gap10_dedx_v_v;
    std::vector<float> _shr_tkfit_gap10_dedx_y_v;

    std::vector<int> _shr_tkfit_dedx_nhits_u_v;
    std::vector<int> _shr_tkfit_dedx_nhits_v_v;
    std::vector<int> _shr_tkfit_dedx_nhits_y_v;

    std::vector<float> _shr_moliere_avg_v;
    std::vector<float> _shr_moliere_rms_v;

    std::vector<float> _shr_llr_pid_u_v;
    std::vector<float> _shr_llr_pid_v_v;
    std::vector<float> _shr_llr_pid_y_v;
    std::vector<float> _shr_llr_pid_v;
    std::vector<float> _shr_llr_pid_score_v;
};

ShowerAnalysis::ShowerAnalysis(const fhicl::ParameterSet &p)
{
    fTRKproducer = p.get<art::InputTag>("TRKproducer", "");
    fCALproducer = p.get<art::InputTag>("CALproducer", "");

    fBacktrackTag = p.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
    fHproducer = p.get<art::InputTag>("Hproducer", "gaushit");
    fEnergyThresholdForMCHits = p.get<float>("EnergyThresholdForMCHits", 0.1);

    fdEdxcmSkip = p.get<float>("dEdxcmSkip", 0.0);
    fdEdxcmLen = p.get<float>("dEdxcmLen", 4.0);
    fLocaldEdx = p.get<bool>("LocaldEdx", false);

    fADCtoE = p.get<std::vector<float>>("ADCtoE");

    fRecalibrateHits = p.get<bool>("RecalibrateHits", false);

    llr_pid_calculator.set_dedx_binning(0, electronphoton_parameters.dedx_edges_pl_0);
    llr_pid_calculator.set_par_binning(0, electronphoton_parameters.parameters_edges_pl_0);
    llr_pid_calculator.set_lookup_tables(0, electronphoton_parameters.dedx_pdf_pl_0);

    llr_pid_calculator.set_dedx_binning(1, electronphoton_parameters.dedx_edges_pl_1);
    llr_pid_calculator.set_par_binning(1, electronphoton_parameters.parameters_edges_pl_1);
    llr_pid_calculator.set_lookup_tables(1, electronphoton_parameters.dedx_pdf_pl_1);

    llr_pid_calculator.set_dedx_binning(2, electronphoton_parameters.dedx_edges_pl_2);
    llr_pid_calculator.set_par_binning(2, electronphoton_parameters.parameters_edges_pl_2);
    llr_pid_calculator.set_lookup_tables(2, electronphoton_parameters.dedx_pdf_pl_2);

    if (fRecalibrateHits)
    {
        llr_pid_calculator.set_corr_par_binning(0, correction_parameters.parameter_correction_edges_pl_0);
        llr_pid_calculator.set_correction_tables(0, correction_parameters.correction_table_pl_0);

        llr_pid_calculator.set_corr_par_binning(1, correction_parameters.parameter_correction_edges_pl_1);
        llr_pid_calculator.set_correction_tables(1, correction_parameters.correction_table_pl_1);

        llr_pid_calculator.set_corr_par_binning(2, correction_parameters.parameter_correction_edges_pl_2);
        llr_pid_calculator.set_correction_tables(2, correction_parameters.correction_table_pl_2);
    }
}

void ShowerAnalysis::configure(fhicl::ParameterSet const &p)
{
}

void ShowerAnalysis::analyseEvent(art::Event const &e, bool fData)
{
    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();
    std::cout << "[ShowerAnalysis::analyseEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: " << _evt << std::endl;
}

void ShowerAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
    art::InputTag clusproducer("pandora");

    common::ProxyCaloColl_t const *tkcalo_proxy = NULL;
    if (fTRKproducer != "")
    {
        tkcalo_proxy = new common::ProxyCaloColl_t(proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer, proxy::withAssociated<anab::Calorimetry>(fCALproducer)));
    }

    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;

    if (!fData)
    {
        art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
        assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
    }

    TVector3 nuvtx;

    for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
    {
        auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

        if (PDG == 12 || PDG == 14)
        {
            double xyz[3] = {};

            auto vtx = slice_pfp_v[i_pfp].get<recob::Vertex>();
            if (vtx.size() != 1)
            {
                std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
            }
            else
            {
                vtx.at(0)->XYZ(xyz);
                nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
            }

            break;
        }
    }

    for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
    {
        auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

        if (PDG == 12 || PDG == 14)
            continue;

        size_t n_shw = slice_pfp_v[i_pfp].get<recob::Shower>().size();
        if (n_shw == 0)
        {
            fillDefault();
        }
        else if (n_shw == 1)
        {
            auto shr = slice_pfp_v[i_pfp].get<recob::Shower>()[0];
            auto const& dedx_v = shr->dEdx();
            auto const& energy_v = shr->Energy();

            if (dedx_v.size() >= 3) {
                _shr_dedx_u_v.push_back(dedx_v[0]);
                _shr_dedx_v_v.push_back(dedx_v[1]);
                _shr_dedx_y_v.push_back(dedx_v[2]);
            } else {
                std::cout << "Warning: Shower dEdx vector has size " << dedx_v.size() 
                          << " for shower ID " << shr->ID() << ", using default values." << std::endl;
                _shr_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
                _shr_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
                _shr_dedx_y_v.push_back(std::numeric_limits<float>::lowest());
            }

            if (energy_v.size() >= 3) {
                _shr_energy_u_v.push_back(energy_v[0]);
                _shr_energy_v_v.push_back(energy_v[1]);
                _shr_energy_y_v.push_back(energy_v[2]);
            } else {
                std::cout << "Warning: Shower energy vector has size " << energy_v.size() 
                          << " for shower ID " << shr->ID() << ", using default values." << std::endl;
                _shr_energy_u_v.push_back(std::numeric_limits<float>::lowest());
                _shr_energy_v_v.push_back(std::numeric_limits<float>::lowest());
                _shr_energy_y_v.push_back(std::numeric_limits<float>::lowest());
            }

            _shr_pfp_id_v.push_back(i_pfp);
            _shr_openangle_v.push_back(shr->OpenAngle());
            _shr_phi_v.push_back(shr->Direction().Phi());
            _shr_theta_v.push_back(shr->Direction().Theta());

            _shr_pitch_u_v.push_back(common::getPitch(shr->Direction().Y(), shr->Direction().Z(), 0));
            _shr_pitch_v_v.push_back(common::getPitch(shr->Direction().Y(), shr->Direction().Z(), 1));
            _shr_pitch_y_v.push_back(common::getPitch(shr->Direction().Y(), shr->Direction().Z(), 2));

            _shr_start_x_v.push_back(shr->ShowerStart().X());
            _shr_start_y_v.push_back(shr->ShowerStart().Y());
            _shr_start_z_v.push_back(shr->ShowerStart().Z());

            _shr_start_U_v.push_back(common::YZtoPlanecoordinate(shr->ShowerStart().Y(), shr->ShowerStart().Z(), 0));
            _shr_start_V_v.push_back(common::YZtoPlanecoordinate(shr->ShowerStart().Y(), shr->ShowerStart().Z(), 1));

            _shr_tkfit_nhits_v.push_back(std::numeric_limits<int>::lowest());
            _shr_tkfit_start_x_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_start_y_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_start_z_v.push_back(std::numeric_limits<float>::lowest());

            _shr_tkfit_start_U_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_start_V_v.push_back(std::numeric_limits<float>::lowest());

            _shr_tkfit_phi_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_theta_v.push_back(std::numeric_limits<float>::lowest());

            _shr_tkfit_pitch_u_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_pitch_v_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_pitch_y_v.push_back(std::numeric_limits<float>::lowest());

            _shr_tkfit_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

            _shr_tkfit_gap10_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_gap10_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
            _shr_tkfit_gap10_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

            _shr_tkfit_dedx_nhits_u_v.push_back(std::numeric_limits<int>::lowest());
            _shr_tkfit_dedx_nhits_v_v.push_back(std::numeric_limits<int>::lowest());
            _shr_tkfit_dedx_nhits_y_v.push_back(std::numeric_limits<int>::lowest());

            _shr_llr_pid_u_v.push_back(0);
            _shr_llr_pid_v_v.push_back(0);
            _shr_llr_pid_y_v.push_back(0);
            _shr_llr_pid_v.push_back(0);
            _shr_llr_pid_score_v.push_back(0);

            _shr_px_v.push_back(shr->Direction().X());
            _shr_py_v.push_back(shr->Direction().Y());
            _shr_pz_v.push_back(shr->Direction().Z());

            _shr_dist_v.push_back((shr->ShowerStart() - nuvtx).Mag());

            float _shrmoliereavg;
            float _shrmoliererms;
            common::GetMoliereRadius(slice_pfp_v[i_pfp], _shrmoliereavg, _shrmoliererms);
            _shr_moliere_rms_v.push_back(_shrmoliererms);
            _shr_moliere_avg_v.push_back(_shrmoliereavg);

            if (tkcalo_proxy == NULL)
                continue;

            for (const common::ProxyCaloElem_t &tk : *tkcalo_proxy)
            {
                if (tk->ID() != int(slice_pfp_v[i_pfp].index()))
                    continue;

                _shr_tkfit_nhits_v.back() = tk->CountValidPoints();
                _shr_tkfit_start_x_v.back() = tk->Start().X();
                _shr_tkfit_start_y_v.back() = tk->Start().Y();
                _shr_tkfit_start_z_v.back() = tk->Start().Z();

                _shr_tkfit_start_U_v.back() = common::YZtoPlanecoordinate(_shr_tkfit_start_y_v.back(), _shr_tkfit_start_z_v.back(), 0);
                _shr_tkfit_start_V_v.back() = common::YZtoPlanecoordinate(_shr_tkfit_start_y_v.back(), _shr_tkfit_start_z_v.back(), 1);

                _shr_tkfit_phi_v.back() = tk->StartDirection().Phi();
                _shr_tkfit_theta_v.back() = tk->StartDirection().Theta();

                _shr_tkfit_pitch_u_v.push_back(common::getPitch(tk->StartDirection().Y(), tk->StartDirection().Z(), 0));
                _shr_tkfit_pitch_v_v.push_back(common::getPitch(tk->StartDirection().Y(), tk->StartDirection().Z(), 1));
                _shr_tkfit_pitch_y_v.push_back(common::getPitch(tk->StartDirection().Y(), tk->StartDirection().Z(), 2));

                auto const tkcalos = tk.get<anab::Calorimetry>();

                float calodEdx;
                int caloNpts;

                for (const auto &tkcalo : tkcalos)
                {
                    auto const& plane = tkcalo->PlaneID().Plane;
                    if (plane > 2)
                        continue;

                    std::vector<float> dqdx_values_corrected;

                    if (fData || !fRecalibrateHits)
                    {
                        if (!fLocaldEdx)
                            dqdx_values_corrected = tkcalo->dQdx();
                        else
                            dqdx_values_corrected = tkcalo->dEdx();
                    }
                    else
                        dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(tkcalo, *tk, assocMCPart, fRecalibrateHits, fEnergyThresholdForMCHits, fLocaldEdx);

                    auto const& xyz_v = tkcalo->XYZ();

                    std::vector<float> x_v, y_v, z_v;
                    std::vector<float> dist_from_start_v;

                    float shr_tkfit_start_sce[3];
                    common::ApplySCECorrectionXYZ(_shr_tkfit_start_x_v.back(),_shr_tkfit_start_y_v.back(),_shr_tkfit_start_z_v.back(),shr_tkfit_start_sce);

                    for (auto xyz : xyz_v)
                    {
                        x_v.push_back(xyz.X());
                        y_v.push_back(xyz.Y());
                        z_v.push_back(xyz.Z());

                        float dist_from_start = common::distance3d(xyz.X(), xyz.Y(), xyz.Z(),
                                                              shr_tkfit_start_sce[0], shr_tkfit_start_sce[1], shr_tkfit_start_sce[2]);
                        dist_from_start_v.push_back(dist_from_start);
                    }

                    std::vector<float> dedx_v;
                    if (!fLocaldEdx)
                    {
                        dedx_v = common::GetdEdxfromdQdx(dqdx_values_corrected, x_v, y_v, z_v, 2.1, fADCtoE[plane]);
                    }
                    else
                    {
                        dedx_v = dqdx_values_corrected;
                    }

                    common::GetTrackFitdEdx(dedx_v, tkcalo->ResidualRange(), fdEdxcmSkip, fdEdxcmLen, calodEdx, caloNpts);

                    if (plane == 0)
                    {
                        _shr_tkfit_dedx_u_v.back() = calodEdx;
                        _shr_tkfit_dedx_nhits_u_v.back() = caloNpts;
                    }

                    if (plane == 1)
                    {
                        _shr_tkfit_dedx_v_v.back() = calodEdx;
                        _shr_tkfit_dedx_nhits_v_v.back() = caloNpts;
                    }

                    if (plane == 2)
                    {
                        _shr_tkfit_dedx_y_v.back() = calodEdx;
                        _shr_tkfit_dedx_nhits_y_v.back() = caloNpts;
                    }
                    common::GetTrackFitdEdx(dedx_v, tkcalo->ResidualRange(), 1.0, fdEdxcmLen, calodEdx, caloNpts);

                    if (plane == 2)
                    {
                        _shr_tkfit_gap10_dedx_y_v.back() = calodEdx;
                    }
                    else if (plane == 1)
                    {
                        _shr_tkfit_gap10_dedx_v_v.back() = calodEdx;
                    }
                    else if (plane == 0)
                    {
                        _shr_tkfit_gap10_dedx_u_v.back() = calodEdx;
                    }

                    std::vector<std::vector<float>> par_values;
                    par_values.push_back(dist_from_start_v);
                    auto const &pitch = tkcalo->TrkPitchVec();
                    par_values.push_back(pitch);

                    float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_v, par_values, plane);
                    if (plane == 0)
                    {
                        _shr_llr_pid_u_v.back() = llr_pid;
                    }
                    else if (plane == 1)
                    {
                        _shr_llr_pid_v_v.back() = llr_pid;
                    }
                    else if (plane == 2)
                    {
                        _shr_llr_pid_y_v.back() = llr_pid;
                    }
                    _shr_llr_pid_v.back() += llr_pid;
                }
                _shr_llr_pid_score_v.back() = atan(_shr_llr_pid_v.back() / 10.) * 2 / 3.14159266;
            }
        }
    }
}

void ShowerAnalysis::setBranches(TTree *_tree)
{
    _tree->Branch("shr_dedx_u_v", "std::vector< float >", &_shr_dedx_u_v);
    _tree->Branch("shr_dedx_v_v", "std::vector< float >", &_shr_dedx_v_v);
    _tree->Branch("shr_dedx_y_v", "std::vector< float >", &_shr_dedx_y_v);

    _tree->Branch("shr_energy_u_v", "std::vector< float >", &_shr_energy_u_v);
    _tree->Branch("shr_energy_v_v", "std::vector< float >", &_shr_energy_v_v);
    _tree->Branch("shr_energy_y_v", "std::vector< float >", &_shr_energy_y_v);

    _tree->Branch("shr_pfp_id_v", "std::vector< size_t >", &_shr_pfp_id_v);

    _tree->Branch("shr_start_x_v", "std::vector< float >", &_shr_start_x_v);
    _tree->Branch("shr_start_y_v", "std::vector< float >", &_shr_start_y_v);
    _tree->Branch("shr_start_z_v", "std::vector< float >", &_shr_start_z_v);
    _tree->Branch("shr_dist_v", "std::vector< float >", &_shr_dist_v);

    _tree->Branch("shr_start_U_v", "std::vector< float >", &_shr_start_U_v);
    _tree->Branch("shr_start_V_v", "std::vector< float >", &_shr_start_V_v);

    _tree->Branch("shr_px_v", "std::vector< float >", &_shr_px_v);
    _tree->Branch("shr_py_v", "std::vector< float >", &_shr_py_v);
    _tree->Branch("shr_pz_v", "std::vector< float >", &_shr_pz_v);

    _tree->Branch("shr_openangle_v", "std::vector< float >", &_shr_openangle_v);
    _tree->Branch("shr_theta_v", "std::vector< float >", &_shr_theta_v);
    _tree->Branch("shr_phi_v", "std::vector< float >", &_shr_phi_v);

    _tree->Branch("shr_pitch_u_v", "std::vector<float>", &_shr_pitch_u_v);
    _tree->Branch("shr_pitch_v_v", "std::vector<float>", &_shr_pitch_v_v);
    _tree->Branch("shr_pitch_y_v", "std::vector<float>", &_shr_pitch_y_v);

    _tree->Branch("shr_tkfit_nhits_v", "std::vector< int >", &_shr_tkfit_nhits_v);
    _tree->Branch("shr_tkfit_start_x_v", "std::vector< float >", &_shr_tkfit_start_x_v);
    _tree->Branch("shr_tkfit_start_y_v", "std::vector< float >", &_shr_tkfit_start_y_v);
    _tree->Branch("shr_tkfit_start_z_v", "std::vector< float >", &_shr_tkfit_start_z_v);

    _tree->Branch("shr_tkfit_start_U_v", "std::vector< float >", &_shr_tkfit_start_U_v);
    _tree->Branch("shr_tkfit_start_V_v", "std::vector< float >", &_shr_tkfit_start_V_v);

    _tree->Branch("shr_tkfit_theta_v", "std::vector< float >", &_shr_tkfit_theta_v);
    _tree->Branch("shr_tkfit_phi_v", "std::vector< float >", &_shr_tkfit_phi_v);

    _tree->Branch("shr_tkfit_pitch_u_v", "std::vector<float>", &_shr_tkfit_pitch_u_v);
    _tree->Branch("shr_tkfit_pitch_v_v", "std::vector<float>", &_shr_tkfit_pitch_v_v);
    _tree->Branch("shr_tkfit_pitch_y_v", "std::vector<float>", &_shr_tkfit_pitch_y_v);

    _tree->Branch("shr_tkfit_dedx_u_v", "std::vector< float >", &_shr_tkfit_dedx_u_v);
    _tree->Branch("shr_tkfit_dedx_v_v", "std::vector< float >", &_shr_tkfit_dedx_v_v);
    _tree->Branch("shr_tkfit_dedx_y_v", "std::vector< float >", &_shr_tkfit_dedx_y_v);

    _tree->Branch("shr_tkfit_gap10_dedx_u_v", "std::vector< float >", &_shr_tkfit_gap10_dedx_u_v);
    _tree->Branch("shr_tkfit_gap10_dedx_v_v", "std::vector< float >", &_shr_tkfit_gap10_dedx_v_v);
    _tree->Branch("shr_tkfit_gap10_dedx_y_v", "std::vector< float >", &_shr_tkfit_gap10_dedx_y_v);

    _tree->Branch("shr_tkfit_dedx_nhits_u_v", "std::vector< int >", &_shr_tkfit_dedx_nhits_u_v);
    _tree->Branch("shr_tkfit_dedx_nhits_v_v", "std::vector< int >", &_shr_tkfit_dedx_nhits_v_v);
    _tree->Branch("shr_tkfit_dedx_nhits_y_v", "std::vector< int >", &_shr_tkfit_dedx_nhits_y_v);

    _tree->Branch("shr_llr_pid_u_v", "std::vector<float>", &_shr_llr_pid_u_v);
    _tree->Branch("shr_llr_pid_v_v", "std::vector<float>", &_shr_llr_pid_v_v);
    _tree->Branch("shr_llr_pid_y_v", "std::vector<float>", &_shr_llr_pid_y_v);
    _tree->Branch("shr_llr_pid_v", "std::vector<float>", &_shr_llr_pid_v);
    _tree->Branch("shr_llr_pid_score_v", "std::vector<float>", &_shr_llr_pid_score_v);

    _tree->Branch("shr_moliere_avg_v", "std::vector< float >", &_shr_moliere_avg_v);
    _tree->Branch("shr_moliere_rms_v", "std::vector< float >", &_shr_moliere_rms_v);
}

void ShowerAnalysis::fillDefault()
{
    _shr_energy_u_v.push_back(std::numeric_limits<float>::lowest());
    _shr_energy_v_v.push_back(std::numeric_limits<float>::lowest());
    _shr_energy_y_v.push_back(std::numeric_limits<float>::lowest());

    _shr_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
    _shr_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
    _shr_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

    _shr_pfp_id_v.push_back(std::numeric_limits<int>::lowest());

    _shr_start_x_v.push_back(std::numeric_limits<float>::lowest());
    _shr_start_y_v.push_back(std::numeric_limits<float>::lowest());
    _shr_start_z_v.push_back(std::numeric_limits<float>::lowest());

    _shr_start_U_v.push_back(std::numeric_limits<float>::lowest());
    _shr_start_V_v.push_back(std::numeric_limits<float>::lowest());

    _shr_openangle_v.push_back(std::numeric_limits<float>::lowest());
    _shr_theta_v.push_back(std::numeric_limits<float>::lowest());
    _shr_phi_v.push_back(std::numeric_limits<float>::lowest());

    _shr_pitch_u_v.push_back(std::numeric_limits<float>::lowest());
    _shr_pitch_v_v.push_back(std::numeric_limits<float>::lowest());
    _shr_pitch_y_v.push_back(std::numeric_limits<float>::lowest());

    _shr_dist_v.push_back(std::numeric_limits<float>::lowest());

    _shr_px_v.push_back(std::numeric_limits<float>::lowest());
    _shr_py_v.push_back(std::numeric_limits<float>::lowest());
    _shr_pz_v.push_back(std::numeric_limits<float>::lowest());

    _shr_tkfit_nhits_v.push_back(std::numeric_limits<int>::lowest());
    _shr_tkfit_start_x_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_start_y_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_start_z_v.push_back(std::numeric_limits<float>::lowest());

    _shr_tkfit_start_U_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_start_V_v.push_back(std::numeric_limits<float>::lowest());

    _shr_tkfit_theta_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_phi_v.push_back(std::numeric_limits<float>::lowest());

    _shr_tkfit_pitch_u_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_pitch_v_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_pitch_y_v.push_back(std::numeric_limits<float>::lowest());

    _shr_tkfit_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

    _shr_tkfit_gap10_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_gap10_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
    _shr_tkfit_gap10_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

    _shr_tkfit_dedx_nhits_u_v.push_back(std::numeric_limits<int>::lowest());
    _shr_tkfit_dedx_nhits_v_v.push_back(std::numeric_limits<int>::lowest());
    _shr_tkfit_dedx_nhits_y_v.push_back(std::numeric_limits<int>::lowest());

    _shr_llr_pid_u_v.push_back(std::numeric_limits<float>::lowest());
    _shr_llr_pid_v_v.push_back(std::numeric_limits<float>::lowest());
    _shr_llr_pid_y_v.push_back(std::numeric_limits<float>::lowest());
    _shr_llr_pid_v.push_back(std::numeric_limits<float>::lowest());
    _shr_llr_pid_score_v.push_back(std::numeric_limits<float>::lowest());

    _shr_moliere_rms_v.push_back(std::numeric_limits<float>::lowest());
    _shr_moliere_avg_v.push_back(std::numeric_limits<float>::lowest());
}

void ShowerAnalysis::resetTTree(TTree *_tree)
{
    _shr_energy_u_v.clear();
    _shr_energy_v_v.clear();
    _shr_energy_y_v.clear();

    _shr_dedx_u_v.clear();
    _shr_dedx_v_v.clear();
    _shr_dedx_y_v.clear();

    _shr_pfp_id_v.clear();

    _shr_start_x_v.clear();
    _shr_start_y_v.clear();
    _shr_start_z_v.clear();

    _shr_start_U_v.clear();
    _shr_start_V_v.clear();

    _shr_openangle_v.clear();
    _shr_theta_v.clear();
    _shr_phi_v.clear();

    _shr_pitch_u_v.clear();
    _shr_pitch_v_v.clear();
    _shr_pitch_y_v.clear();

    _shr_dist_v.clear();

    _shr_px_v.clear();
    _shr_py_v.clear();
    _shr_pz_v.clear();

    _shr_tkfit_nhits_v.clear();
    _shr_tkfit_start_x_v.clear();
    _shr_tkfit_start_y_v.clear();
    _shr_tkfit_start_z_v.clear();

    _shr_tkfit_start_U_v.clear();
    _shr_tkfit_start_V_v.clear();

    _shr_tkfit_theta_v.clear();
    _shr_tkfit_phi_v.clear();

    _shr_tkfit_pitch_u_v.clear();
    _shr_tkfit_pitch_v_v.clear();
    _shr_tkfit_pitch_y_v.clear();

    _shr_tkfit_dedx_u_v.clear();
    _shr_tkfit_dedx_v_v.clear();
    _shr_tkfit_dedx_y_v.clear();

    _shr_tkfit_gap10_dedx_u_v.clear();
    _shr_tkfit_gap10_dedx_v_v.clear();
    _shr_tkfit_gap10_dedx_y_v.clear();

    _shr_tkfit_dedx_nhits_u_v.clear();
    _shr_tkfit_dedx_nhits_v_v.clear();
    _shr_tkfit_dedx_nhits_y_v.clear();

    _shr_llr_pid_u_v.clear();
    _shr_llr_pid_v_v.clear();
    _shr_llr_pid_y_v.clear();
    _shr_llr_pid_v.clear();
    _shr_llr_pid_score_v.clear();

    _shr_moliere_rms_v.clear();
    _shr_moliere_avg_v.clear();
}

DEFINE_ART_CLASS_TOOL(ShowerAnalysis)
}

#endif