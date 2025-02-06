#ifndef EVENTCLASSIFIER_H
#define EVENTCLASSIFIER_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"

#include <vector>
#include <string>
#include <limits>

namespace signature {

enum EventType {
    kSignal = 0,
    kBeamNeutrino,
    kCosmicRay,
    kExternal,
    kOther
};

class EventClassifier {
public:
    EventClassifier(fhicl::ParameterSet const& pset);

    EventType classifyEvent(const art::Event& e);
    const Pattern& getPattern(const art::Event& e);

private:
    simb::Origin_t getTruthOrigin(const art::Event& e);
    void createPattern(const art::Event& e);
    bool isContained(const double point[3]) const;
    bool isSignal(const art::Event& e);

    float _min_muon_momentum;
    float _min_pion_momentum;
    std::vector<float> _fiducial_offsets; // {x_start, y_start, z_start, x_end, y_end, z_end}

    art::InputTag _MCTproducer;

    std::vector<std::unique_ptr<SignatureToolBase>> _signatureToolsVec;
    Pattern _pattern;
    bool _pattern_found;
};

EventClassifier::EventClassifier(fhicl::ParameterSet const& pset)
    : _min_muon_momentum{pset.get<float>("MinMuonMomentum", 0.0)},
      _min_pion_momentum{pset.get<float>("MinPionMomentum", 0.0)},
      _fiducial_offsets{pset.get<std::vector<float>>("FiducialOffsets", {10, 10, 10, 10, 10, 10})},
      _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
{
    const auto& tools_pset = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (const auto& tool_label : tools_pset.get_pset_names()) {
        auto const& tool_pset = tools_pset.get<fhicl::ParameterSet>(tool_label);
        _signatureToolsVec.push_back(art::make_tool<SignatureToolBase>(tool_pset));
    }
}

simb::Origin_t EventClassifier::getTruthOrigin(const art::Event& e) {
    art::Handle<std::vector<simb::MCTruth>> truthHandle;
    e.getByLabel(_MCTproducer, truthHandle);

    if (!truthHandle.isValid())
        return simb::kUnknown;

    return truthHandle->front().Origin();
}

void EventClassifier::createPattern(const art::Event& e) {
    _pattern.clear();  
    _pattern_found = true;

    for (auto& tool : _signatureToolsVec) {
        Signature signature;
        if (tool->constructSignature(e, signature)) {
            _pattern.emplace_back(tool->getSignatureType(), signature);  
        } else {
            _pattern_found = false;
        }
    }
}

const Pattern& EventClassifier::getPattern(const art::Event& e) {
    if (!_pattern_found) 
        this->createPattern(e);

    return _pattern;
}

bool EventClassifier::isContained(const double point[3]) const {
    art::ServiceHandle<geo::Geometry> geo;
    geo::TPCGeo const &tpc = geo->TPC();
    geo::BoxBoundedGeo bound_box = tpc.ActiveBoundingBox();

    std::vector<double> bounds = {bound_box.MinX(), bound_box.MaxX(),
                                  bound_box.MinY(), bound_box.MaxY(),
                                  bound_box.MinZ(), bound_box.MaxZ()};

    bool is_x = point[0] > (bounds[0] + _fiducial_offsets[0]) && point[0] < (bounds[1] - _fiducial_offsets[3]);
    bool is_y = point[1] > (bounds[2] + _fiducial_offsets[1]) && point[1] < (bounds[3] - _fiducial_offsets[4]);
    bool is_z = point[2] > (bounds[4] + _fiducial_offsets[2]) && point[2] < (bounds[5] - _fiducial_offsets[5]);

    return is_x && is_y && is_z;
}

bool EventClassifier::isSignal(const art::Event& e) {
    art::Handle<std::vector<simb::MCTruth>> truthHandle;
    e.getByLabel(_MCTproducer, truthHandle);

    if (!truthHandle.isValid()) 
        return false;

    this->createPattern(e);
    if (_pattern_found) 
        return false; 

    const simb::MCTruth& truth = truthHandle->front();
    const auto& neutrino = truth.GetNeutrino();

    const double point[3] = {neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz()};
    if (!this->isContained(point)) {
        return false;
    }

    for (const auto& [_, signature] : _pattern) {
        for (const auto& particle : signature) {
            int pdg_code = particle->PdgCode();
            float momentum = particle->P();
            std::string process = particle->Process();

            if (std::abs(pdg_code) == 13 && process == "primary" && momentum < _min_muon_momentum) 
                return false;  // Muon
            if (std::abs(pdg_code) == 211 && process == "Decay" && momentum < _min_pion_momentum) 
                return false;  // Pions
        }
    }

    return true;  
}

EventType EventClassifier::classifyEvent(const art::Event& e) {
    simb::Origin_t origin = this->getTruthOrigin(e);

    if (origin == simb::kBeamNeutrino && this->isSignal(e)) 
        return kSignal;

    switch (origin) {
        case simb::kBeamNeutrino:       return kBeamNeutrino;
        case simb::kCosmicRay:          return kCosmicRay;
        case simb::kSingleParticle:     return kExternal;
        default:                        return kOther;
    }
}

} 

#endif // EVENTCLASSIFIER_H
