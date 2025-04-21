#ifndef SIGNALDEFINITION_ANALYSIS_CXX
#define SIGNALDEFINITION_ANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "TTree.h"
#include "TVector3.h"
#include <memory>
#include <vector>
#include <array>
#include <optional>
#include <unordered_set>
#include <map>
#include <algorithm>

namespace analysis 
{
    enum SignatureType {
        kEmptySignature = 0,
        kCosmicSignature,
        kNeutrinoSignature,
        kPrimaryMuonSignature,
        kChargedKaonSignature,
        kKaonShortSignature,
        kLambdaSignature,
        kChargedSigmaSignature
    };

    enum class EventType {
        kSignal = 0,
        kBeamNeutrino,
        kCosmicRay,
        kOther
    };

    using Signature = std::vector<art::Ptr<simb::MCParticle>>;
    using Pattern = std::vector<std::pair<SignatureType, Signature>>;
    using ParticleMap = std::map<int, art::Ptr<simb::MCParticle>>;

    class SignalDefinitionAnalysis : public AnalysisToolBase {
    public:
        explicit SignalDefinitionAnalysis(const fhicl::ParameterSet& pset);
        ~SignalDefinitionAnalysis() = default;

        SignalDefinitionAnalysis(const SignalDefinitionAnalysis&) = delete;
        SignalDefinitionAnalysis(SignalDefinitionAnalysis&&) = delete;
        SignalDefinitionAnalysis& operator=(const SignalDefinitionAnalysis&) = delete;
        SignalDefinitionAnalysis& operator=(SignalDefinitionAnalysis&&) = delete;

        void configure(const fhicl::ParameterSet& pset) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;
        void analyseEvent(const art::Event& event, bool is_data) override;
        void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slicePfpVector, bool is_data, bool selected) override {}

    private:
        art::InputTag _MCTproducer;
        art::InputTag _MCRproducer;
        art::InputTag _MCPproducer;
        art::InputTag _HitProducer;
        art::InputTag _BacktrackTag;

        std::vector<float> _fiducialOffsets;
        double _hit_exclusivity_thresh;
        double _sig_exclusivity_thresh;

        art::ServiceHandle<geo::Geometry> _geoService;

        int event_type_;
        int is_signal_;
        std::vector<int> signature_types_;
        std::vector<int> clarity_U_, clarity_V_, clarity_W_;
        int pass_clarity_;
        std::vector<int> particle_pdg_;
        std::vector<double> particle_energy_;
        std::vector<double> particle_px_, particle_py_, particle_pz_;
        std::array<std::vector<double>, static_cast<size_t>(PandoraView::N_VIEWS)> exclusivity_ratio_;
        std::array<std::vector<double>, static_cast<size_t>(PandoraView::N_VIEWS)> hit_exclusivity_fraction_;
        std::vector<double> vertex_distance_;
        float _mcf_nu_e;
        float _mcf_lep_e;
        int _mcf_actvol;
        int _mcf_nmm, _mcf_nmp, _mcf_nem, _mcf_nep, _mcf_np0, _mcf_npp, _mcf_npm;
        int _mcf_nkp, _mcf_nkm, _mcf_nk0, _mcf_npr, _mcf_nne;
        float _mcf_mcshr_elec_etot;
        int _mcf_pass_ccpi0, _mcf_pass_ncpi0, _mcf_pass_ccnopi, _mcf_pass_ncnopi;
        int _mcf_pass_cccpi, _mcf_pass_nccpi;

        std::optional<TVector3> getPrimaryVertex(const art::Event& e) {
            art::Handle<std::vector<simb::MCTruth>> MCThandle;
            e.getByLabel(_MCTproducer, MCThandle);
            if (MCThandle.isValid() && !MCThandle->empty()) {
                const auto& mct = MCThandle->at(0);
                const auto& nu = mct.GetNeutrino();
                return TVector3(nu.Nu().Vx(), nu.Nu().Vy(), nu.Nu().Vz());
            }
            return std::nullopt;
        }

        bool isContained(const art::Event& e) {
            art::Handle<std::vector<simb::MCTruth>> truthHandle;
            e.getByLabel(_MCTproducer, truthHandle);
            if (!truthHandle.isValid() || truthHandle->empty()) return false;
            const simb::MCTruth& truth = truthHandle->front();
            const auto& neutrino = truth.GetNeutrino();
            double point[3] = {neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz()};
            geo::TPCGeo const& tpc = _geoService->TPC();
            geo::BoxBoundedGeo boundBox = tpc.ActiveBoundingBox();
            std::vector<double> bounds = {boundBox.MinX(), boundBox.MaxX(), boundBox.MinY(), boundBox.MaxY(), boundBox.MinZ(), boundBox.MaxZ()};
            bool isX = point[0] > (bounds[0] + _fiducialOffsets[0]) && point[0] < (bounds[1] - _fiducialOffsets[3]);
            bool isY = point[1] > (bounds[2] + _fiducialOffsets[1]) && point[1] < (bounds[3] - _fiducialOffsets[4]);
            bool isZ = point[2] > (bounds[4] + _fiducialOffsets[2]) && point[2] < (bounds[5] - _fiducialOffsets[5]);
            return isX && isY && isZ;
        }

        simb::Origin_t getTruthOrigin(const art::Event& e) {
            art::Handle<std::vector<simb::MCTruth>> truthHandle;
            e.getByLabel(_MCTproducer, truthHandle);
            if (!truthHandle.isValid() || truthHandle->empty()) return simb::kUnknown;
            return truthHandle->front().Origin();
        }

        std::vector<art::Ptr<simb::MCParticle>> getDaughters(const art::Ptr<simb::MCParticle>& particle, const std::map<int, art::Ptr<simb::MCParticle>>& mcp_map) {
            std::vector<art::Ptr<simb::MCParticle>> daughters;
            for (int i = 0; i < particle->NumberDaughters(); ++i) {
                int daughter_id = particle->Daughter(i);
                if (mcp_map.count(daughter_id)) {
                    daughters.push_back(mcp_map.at(daughter_id));
                }
            }
            return daughters;
        }

        bool matchesDecayProducts(const art::Ptr<simb::MCParticle>& particle, const std::map<int, art::Ptr<simb::MCParticle>>& mcp_map, const std::vector<int>& pdg_codes) {
            auto daughters = getDaughters(particle, mcp_map);
            if (daughters.size() != pdg_codes.size()) return false;
            std::vector<int> daughter_pdgs;
            for (const auto& daughter : daughters) {
                daughter_pdgs.push_back(daughter->PdgCode());
            }
            std::sort(daughter_pdgs.begin(), daughter_pdgs.end());
            std::vector<int> expected_pdgs = pdg_codes;
            std::sort(expected_pdgs.begin(), expected_pdgs.end());
            return daughter_pdgs == expected_pdgs;
        }

        void findKaonShortSignature(const art::Event& e, Signature& signature, bool& found) {
            auto mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
            std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
            for (size_t i = 0; i < mcp_h->size(); ++i) {
                mcp_map[mcp_h->at(i).TrackId()] = art::Ptr<simb::MCParticle>(mcp_h, i);
            }

            for (const auto& mcp : *mcp_h) {
                if (abs(mcp.PdgCode()) == 311 && mcp.Process() == "primary" && mcp.EndProcess() == "Decay" && !found) {
                    auto neutral_kaon_daughters = getDaughters(mcp_map.at(mcp.TrackId()), mcp_map);
                    for (const auto& kaon_shorts : neutral_kaon_daughters) {
                        if (kaon_shorts->PdgCode() == 310 && kaon_shorts->EndProcess() == "Decay") {
                            if (matchesDecayProducts(kaon_shorts, mcp_map, {211, -211})) {
                                auto decay_products = getDaughters(kaon_shorts, mcp_map);
                                found = true;
                                for (const auto& pion : decay_products) {
                                    signature.push_back(pion);
                                    auto addDaughterInteractions = [&](const art::Ptr<simb::MCParticle>& particle, auto& self) -> void {
                                        auto daughters = getDaughters(particle, mcp_map);
                                        for (const auto& daugh : daughters) {
                                            if (daugh->PdgCode() == particle->PdgCode()) {
                                                signature.push_back(daugh);
                                                self(daugh, self);
                                            }
                                        }
                                    };
                                    addDaughterInteractions(pion, addDaughterInteractions);
                                }
                                if (found) return;
                            }
                        }
                    }
                }
            }
        }

        TVector3 findKaonShortVertex(const art::Event& e, const Signature& signature) {
            if (!signature.empty()) {
                const auto& first_particle = signature[0];
                return TVector3(first_particle->Position(0).X(), first_particle->Position(0).Y(), first_particle->Position(0).Z());
            }
            return TVector3();
        }

        bool evaluateHitExclusivity(const art::Event& e, const Signature& sig, PandoraView view, double& exclusivity_ratio, double& hit_exclusivity_fraction) {
            art::Handle<std::vector<recob::Hit>> hit_handle;
            e.getByLabel(_HitProducer, hit_handle);
            if (!hit_handle.isValid()) return false;

            art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcp_bkth_assoc(hit_handle, e, _BacktrackTag);
            if (!mcp_bkth_assoc.isValid()) return false;

            std::unordered_set<int> sig_tracks;
            for (const auto& particle : sig) {
                sig_tracks.insert(particle->TrackId());
            }

            double exclusive_q = 0.0;
            double inclusive_q = 0.0;
            std::unordered_set<size_t> ass_hits;
            std::unordered_set<size_t> exc_hits;

            for (size_t i = 0; i < hit_handle->size(); ++i) {
                const recob::Hit& hit = (*hit_handle)[i];
                if (hit.View() != static_cast<int>(view)) continue;

                size_t hit_key = i;
                const auto& ass_part = mcp_bkth_assoc.at(hit_key);
                const auto& ass_data = mcp_bkth_assoc.data(hit_key);
                bool is_assoc = false;
                bool is_excus = false;

                for (size_t j = 0; j < ass_part.size(); ++j) {
                    if (sig_tracks.count(ass_part[j]->TrackId())) {
                        is_assoc = true;
                        double charge = ass_data[j]->numElectrons * ass_data[j]->ideNFraction;
                        inclusive_q += charge;
                        if (ass_data[j]->ideNFraction > _hit_exclusivity_thresh) {
                            exclusive_q += charge;
                            is_excus = true;
                            break;
                        }
                    }
                }
                if (is_assoc) {
                    ass_hits.insert(hit_key);
                    if (is_excus) {
                        exc_hits.insert(hit_key);
                    }
                }
            }

            exclusivity_ratio = inclusive_q > 0 ? exclusive_q / inclusive_q : 0.0;
            hit_exclusivity_fraction = ass_hits.empty() ? 0.0 : static_cast<double>(exc_hits.size()) / ass_hits.size();

            return exclusivity_ratio > _sig_exclusivity_thresh;
        }
    };

    inline SignalDefinitionAnalysis::SignalDefinitionAnalysis(const fhicl::ParameterSet& pset) {
        configure(pset);
    }

    inline void SignalDefinitionAnalysis::configure(const fhicl::ParameterSet& pset) {
        _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");
        _MCRproducer = pset.get<art::InputTag>("MCRproducer", "largeant");
        _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
        _HitProducer = pset.get<art::InputTag>("HitProducer", "gaushit");
        _BacktrackTag = pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
        _fiducialOffsets = pset.get<std::vector<float>>("FiducialOffsets", {10, 10, 10, 10, 10, 10});
        _hit_exclusivity_thresh = pset.get<double>("HitExclusivityThreshold", 0.50);
        _sig_exclusivity_thresh = pset.get<double>("SignatureExclusivityThreshold", 0.70);
    }

    inline void SignalDefinitionAnalysis::setBranches(TTree* tree) {
        tree->Branch("event_type", &event_type_);
        tree->Branch("is_signal", &is_signal_);
        tree->Branch("signature_types", &signature_types_);
        tree->Branch("clarity_U", &clarity_U_);
        tree->Branch("clarity_V", &clarity_V_);
        tree->Branch("clarity_W", &clarity_W_);
        tree->Branch("pass_clarity", &pass_clarity_);
        tree->Branch("exclusivity_ratio_U", &exclusivity_ratio_[0]);
        tree->Branch("exclusivity_ratio_V", &exclusivity_ratio_[1]);
        tree->Branch("exclusivity_ratio_W", &exclusivity_ratio_[2]);
        tree->Branch("hit_exclusivity_fraction_U", &hit_exclusivity_fraction_[0]);
        tree->Branch("hit_exclusivity_fraction_V", &hit_exclusivity_fraction_[1]);
        tree->Branch("hit_exclusivity_fraction_W", &hit_exclusivity_fraction_[2]);
        tree->Branch("particle_pdg", &particle_pdg_);
        tree->Branch("particle_energy", &particle_energy_);
        tree->Branch("particle_px", &particle_px_);
        tree->Branch("particle_py", &particle_py_);
        tree->Branch("particle_pz", &particle_pz_);
        tree->Branch("vertex_distance", &vertex_distance_);
        tree->Branch("mcf_nu_e", &_mcf_nu_e, "mcf_nu_e/F");
        tree->Branch("mcf_lep_e", &_mcf_lep_e, "mcf_lep_e/F");
        tree->Branch("mcf_actvol", &_mcf_actvol, "mcf_actvol/I");
        tree->Branch("mcf_nmm", &_mcf_nmm, "mcf_nmm/I");
        tree->Branch("mcf_nmp", &_mcf_nmp, "mcf_nmp/I");
        tree->Branch("mcf_nem", &_mcf_nem, "mcf_nem/I");
        tree->Branch("mcf_nep", &_mcf_nep, "mcf_nep/I");
        tree->Branch("mcf_np0", &_mcf_np0, "mcf_np0/I");
        tree->Branch("mcf_npp", &_mcf_npp, "mcf_npp/I");
        tree->Branch("mcf_npm", &_mcf_npm, "mcf_npm/I");
        tree->Branch("mcf_nkp", &_mcf_nkp, "mcf_nkp/I");
        tree->Branch("mcf_nkm", &_mcf_nkm, "mcf_nkm/I");
        tree->Branch("mcf_nk0", &_mcf_nk0, "mcf_nk0/I");
        tree->Branch("mcf_npr", &_mcf_npr, "mcf_npr/I");
        tree->Branch("mcf_nne", &_mcf_nne, "mcf_nne/I");
        tree->Branch("mcf_mcshr_elec_etot", &_mcf_mcshr_elec_etot, "mcf_mcshr_elec_etot/F");
        tree->Branch("mcf_pass_ccpi0", &_mcf_pass_ccpi0, "mcf_pass_ccpi0/I");
        tree->Branch("mcf_pass_ncpi0", &_mcf_pass_ncpi0, "mcf_pass_ncpi0/I");
        tree->Branch("mcf_pass_ccnopi", &_mcf_pass_ccnopi, "mcf_pass_ccnopi/I");
        tree->Branch("mcf_pass_ncnopi", &_mcf_pass_ncnopi, "mcf_pass_ncnopi/I");
        tree->Branch("mcf_pass_cccpi", &_mcf_pass_cccpi, "mcf_pass_cccpi/I");
        tree->Branch("mcf_pass_nccpi", &_mcf_pass_nccpi, "mcf_pass_nccpi/I");
    }

    inline void SignalDefinitionAnalysis::resetTTree(TTree* tree) {
        event_type_ = -1;
        is_signal_ = 0;
        signature_types_.clear();
        clarity_U_.clear();
        clarity_V_.clear();
        clarity_W_.clear();
        pass_clarity_ = 0;
        particle_pdg_.clear();
        particle_energy_.clear();
        particle_px_.clear();
        particle_py_.clear();
        particle_pz_.clear();
        for (auto& vec : exclusivity_ratio_) vec.clear();
        for (auto& vec : hit_exclusivity_fraction_) vec.clear();
        vertex_distance_.clear();
        _mcf_nu_e = -1.0;
        _mcf_lep_e = -1.0;
        _mcf_actvol = -1;
        _mcf_nmm = _mcf_nmp = _mcf_nem = _mcf_nep = _mcf_np0 = 0;
        _mcf_npp = _mcf_npm = _mcf_nkp = _mcf_nkm = _mcf_nk0 = 0;
        _mcf_npr = _mcf_nne = 0;
        _mcf_mcshr_elec_etot = -1.0;
        _mcf_pass_ccpi0 = _mcf_pass_ncpi0 = _mcf_pass_ccnopi = -1;
        _mcf_pass_ncnopi = _mcf_pass_cccpi = _mcf_pass_nccpi = -1;
    }

    inline void SignalDefinitionAnalysis::analyseEvent(const art::Event& e, bool is_data) {
        if (is_data) return;

        // Truth information
        auto mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
        auto mct = mct_h->at(0);
        auto neutrino = mct.GetNeutrino();
        auto nu = neutrino.Nu();
        _mcf_nu_e = nu.Trajectory().E(0);
        _mcf_lep_e = neutrino.Lepton().E();

        geo::TPCGeo const& tpc = _geoService->TPC();
        geo::BoxBoundedGeo activeBox = tpc.ActiveBoundingBox();
        double vx = nu.Vx(), vy = nu.Vy(), vz = nu.Vz();
        _mcf_actvol = (activeBox.ContainsPosition(geo::Point_t(vx, vy, vz))) ? 1 : 0;

        size_t npart = mct.NParticles();
        _mcf_nmm = _mcf_nmp = _mcf_nem = _mcf_nep = _mcf_np0 = 0;
        _mcf_npp = _mcf_npm = _mcf_nkp = _mcf_nkm = _mcf_nk0 = 0;
        _mcf_npr = _mcf_nne = 0;
        for (size_t i = 0; i < npart; ++i) {
            auto const& part = mct.GetParticle(i);
            if (part.StatusCode() != 1) continue;
            int pdg = part.PdgCode();
            if (pdg == 13) _mcf_nmm++;
            else if (pdg == -13) _mcf_nmp++;
            else if (pdg == 11) _mcf_nem++;
            else if (pdg == -11) _mcf_nep++;
            else if (pdg == 111) _mcf_np0++;
            else if (pdg == 211) _mcf_npp++;
            else if (pdg == -211) _mcf_npm++;
            else if (pdg == 321) _mcf_nkp++;
            else if (pdg == -321) _mcf_nkm++;
            else if (pdg == 311) _mcf_nk0++;
            else if (pdg == 2212) _mcf_npr++;
            else if (pdg == 2112) _mcf_nne++;
        }

        float maxElecMCShwEMeV = 0.0;
        auto mcshr_h = e.getValidHandle<std::vector<sim::MCShower>>(_MCRproducer);
        for (const auto& mcs : *mcshr_h) {
            if (std::abs(mcs.PdgCode()) == 11) {
                if (mcs.Start().E() > maxElecMCShwEMeV) {
                    _mcf_mcshr_elec_etot = mcs.Start().E();
                    maxElecMCShwEMeV = _mcf_mcshr_elec_etot;
                }
            }
        }

        _mcf_pass_ccpi0 = (_mcf_actvol == 1 && _mcf_nmm == 1 && _mcf_nem == 0 && _mcf_nep == 0 && _mcf_np0 == 1) ? 1 : 0;
        _mcf_pass_ncpi0 = (_mcf_actvol == 1 && _mcf_nmm == 0 && _mcf_nmp == 0 && _mcf_nem == 0 && _mcf_nep == 0 && _mcf_np0 == 1) ? 1 : 0;
        _mcf_pass_ccnopi = (_mcf_actvol == 1 && _mcf_nmm == 1 && _mcf_nem == 0 && _mcf_nep == 0 &&
                            _mcf_np0 == 0 && _mcf_npp == 0 && _mcf_npm == 0 &&
                            (((_mcf_lep_e - 0.105) > 0.02 && _mcf_lep_e < 0.3) || _mcf_mcshr_elec_etot > 15)) ? 1 : 0;
        _mcf_pass_ncnopi = (_mcf_actvol == 1 && _mcf_nmm == 0 && _mcf_nmp == 0 && _mcf_nem == 0 && _mcf_nep == 0 &&
                            _mcf_np0 == 0 && _mcf_npp == 0 && _mcf_npm == 0 && _mcf_nu_e > 0.9) ? 1 : 0;
        _mcf_pass_cccpi = (_mcf_actvol == 1 && _mcf_nmm == 1 && _mcf_nem == 0 && _mcf_nep == 0 &&
                           _mcf_np0 == 0 && (_mcf_npp == 1 || _mcf_npm == 1) &&
                           (((_mcf_lep_e - 0.105) > 0.02 && _mcf_lep_e < 0.4) || _mcf_mcshr_elec_etot > 35)) ? 1 : 0;
        _mcf_pass_nccpi = (_mcf_actvol == 1 && _mcf_nmm == 0 && _mcf_nmp == 0 && _mcf_nem == 0 && _mcf_nep == 0 &&
                           _mcf_np0 == 0 && (_mcf_npp == 1 || _mcf_npm == 1)) ? 1 : 0;

        // Event classification
        simb