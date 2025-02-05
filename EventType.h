#ifndef EVENTTYPES_H
#define EVENTTYPES_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/LArTrainingData.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <iostream>
#include <limits>
#include <string>

class EventTypes {
public:
    enum Topology {
        kTopNue,
        kTopNumu,
        kTopNutauE,
        kTopNutauMu,
        kTopNutauHad,
        kTopNC,
        kTopNueLike,
        kTopNumuLike,
        kTopNutauLike,
        kTopNCLike,
        kUnknown,
        kOther,
        kOOFV
    };

    EventTypes();

    InteractionType GetInteractionType(const simb::MCNeutrino& truth) const;
    InteractionType GetInteractionTypeFromSlice(int pdg, bool iscc, int trueMode) const;
    void GetTopology(const art::Ptr<simb::MCTruth> truth, unsigned int nTopologyHits = 0);


    Topology GetTopologyType() const;
  
    Topology CategorizeEvent(const simb::MCNeutrino& truth, bool isMC, bool vertexInFV, int interactionType, int pdgCode) const;



private:
    void CountFinalStateParticles(const art::Ptr<simb::MCTruth> truth, unsigned int nTopologyHits);

    int nProton, nPion, nPizero, nNeutron, pdgCode, tauMode;
};

EventTypes::EventTypes() 
    : nProton(0), nPion(0), nPizero(0), nNeutron(0), pdgCode(0), tauMode(0) {}

InteractionType EventTypes::GetInteractionType(const simb::MCNeutrino& truth) const {
    int pdg = truth.Nu().PdgCode();
    bool iscc = truth.CCNC() == simb::kCC;
    int trueMode = truth.Mode();

    if (iscc) {
        if (abs(pdg) == 14) {
            switch (trueMode) {
                case simb::kQE: return kNumuQE;
                case simb::kRes: return kNumuRes;
                case simb::kDIS: return kNumuDIS;
                default: return kNumuOther;
            }
        } else if (abs(pdg) == 12) {
            switch (trueMode) {
                case simb::kQE: return kNueQE;
                case simb::kRes: return kNueRes;
                case simb::kDIS: return kNueDIS;
                default: return kNueOther;
            }
        } else if (abs(pdg) == 16) {
            switch (trueMode) {
                case simb::kQE: return kNutauQE;
                case simb::kRes: return kNutauRes;
                case simb::kDIS: return kNutauDIS;
                default: return kNutauOther;
            }
        }
    } else if (trueMode == simb::kNuElectronElastic) {
        return kNuElectronElastic;
    }
    return kNC;
}

void EventTypes::GetTopology(const art::Ptr<simb::MCTruth> truth, unsigned int nTopologyHits) {
    const simb::MCNeutrino& nu = truth->GetNeutrino();

    pdgCode = (nu.CCNC() == simb::kCC) ? nu.Nu().PdgCode() : 1;
    tauMode = kNotNutau;

    if (abs(pdgCode) == 16) {
        tauMode = kNutauHad;
        for (int p = 0; p < truth->NParticles(); ++p) {
            const auto& particle = truth->GetParticle(p);
            if (particle.StatusCode() != 1) continue;

            int pdg = abs(particle.PdgCode());
            int parent = particle.Mother();
            while (parent > 0)
                parent = truth->GetParticle(parent).Mother();

            if (parent == 0) {
                if (pdg == 11) { tauMode = kNutauE; break; }
                else if (pdg == 13) { tauMode = kNutauMu; break; }
            }
        }
    }

    CountFinalStateParticles(truth, nTopologyHits);
}

EventTypes::Topology EventTypes::CategorizeEvent(const simb::MCNeutrino& truth, bool isMC, bool vertexInFV, int interactionType, int pdgCode) const {
    int abs_mc_nu_pdg = std::abs(pdgCode);
    if (!isMC || (abs_mc_nu_pdg != 12 && abs_mc_nu_pdg != 14 && abs_mc_nu_pdg != 16)) {
        return kUnknown;
    }

    if (!vertexInFV) {
        return kOOFV;
    }

    if (truth.CCNC() == simb::kNC) {
        return kTopNC;
    }

    if (pdgCode == 12) {
        return kTopNue;
    } 
    
    if (pdgCode != 14) {
        return kOther;
    }

    bool Is_CC1muNp0pi_Event = (nProton >= 1) && (nPion == 0) && (nPizero == 0) && (nNeutron >= 0);

    if (Is_CC1muNp0pi_Event) {
        if (nProton == 1) {
            if (interactionType == simb::kQE) return kTopNumu;
            else if (interactionType == simb::kRes) return kTopNumuLike;
            else return kOther;
        } else if (nProton == 2) {
            if (interactionType == simb::kQE) return kTopNumu;
            else if (interactionType == simb::kRes) return kTopNumuLike;
            else return kOther;
        } else {
            if (interactionType == simb::kQE) return kTopNumu;
            else if (interactionType == simb::kRes) return kTopNumuLike;
            else return kOther;
        }
    } else if (nPion > 0) {
        return kOther;
    } else if (nProton == 0) {
        if (interactionType == simb::kQE) return kTopNumuLike;
        else if (interactionType == simb::kRes) return kTopNumuLike;
        else return kOther;
    }
    return kOther;
}

void EventTypes::CountFinalStateParticles(const art::Ptr<simb::MCTruth> truth, unsigned int nTopologyHits) {
    nProton = 0;
    nPion = 0;
    nPizero = 0;
    nNeutron = 0;

    art::ServiceHandle<cheat::BackTrackerService> backTrack;
    art::ServiceHandle<cheat::ParticleInventoryService> partService;

    for (auto const& thisPart : partService->MCTruthToParticles_Ps(truth)) {
        const simb::MCParticle& part = *thisPart;

        if (part.StatusCode() != 1 || part.Mother() != 0 || part.PdgCode() > 1000000) continue;

        unsigned int nSimIDE = backTrack->TrackIdToSimIDEs_Ps(part.TrackId()).size();

        if ((part.PdgCode() == 111 || part.PdgCode() == 2112)) {
            for (int d = 0; d < part.NumberDaughters(); ++d) {
                nSimIDE += backTrack->TrackIdToSimIDEs_Ps(part.Daughter(d)).size();
            }
        }

        if (nSimIDE < nTopologyHits) continue;

        switch (abs(part.PdgCode())) {
            case 111: ++nPizero; break;
            case 211: ++nPion; break;
            case 2112: ++nNeutron; break;
            case 2212: ++nProton; break;
            default: break;
        }
    }

}


EventTypes::Topology EventTypes::GetTopologyType() const {
    if (abs(pdgCode) == 12) return kTopNue;
    if (abs(pdgCode) == 14) return kTopNumu;
    if (abs(pdgCode) == 16) {
        if (tauMode == kNutauE) return kTopNutauE;
        if (tauMode == kNutauMu) return kTopNutauMu;
        if (tauMode == kNutauHad) return kTopNutauHad;
    }
    if (pdgCode == 1) return kTopNC;
    throw std::runtime_error("Topology type not recognised!");
}

#endif // EVENTTYPES_H
