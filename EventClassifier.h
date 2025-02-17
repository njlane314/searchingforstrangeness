#ifndef EVENTCLASSIFIER_H
#define EVENTCLASSIFIER_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"

#include "ClarityTools/ClarityToolBase.h"

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
    
    bool passesPlaneClarity(const Signature& sig, common::PandoraView view) const;
    bool passesThreePlaneClarity(const Signature& sig) const;
    bool passClarityFilter() const;

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
    std::vector<std::unique_ptr<ClarityToolBase>> _clarityToolsVec;
    Pattern _pattern;
    bool _pattern_found;

    struct ClarityResult {
        Type type;
        Signature signature;
        std::array<bool, common::N_VIEWS> passes;

        ClarityResult(Type t, const Signature& s)
        : type(t), signature(s)
        {
            passes.fill(true);
        }
    };

    std::vector<ClarityResult> _clarity_results;
};

EventClassifier::EventClassifier(fhicl::ParameterSet const& pset)
    : _min_muon_momentum{pset.get<float>("MinMuonMomentum", 0.0)},
      _min_pion_momentum{pset.get<float>("MinPionMomentum", 0.0)},
      _fiducial_offsets{pset.get<std::vector<float>>("FiducialOffsets", {10, 10, 10, 10, 10, 10})},
      _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
{
    const auto& sig_tools_pset = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (const auto& tool_label : sig_tools_pset.get_pset_names()) {
        auto const& tool_pset = sig_tools_pset.get<fhicl::ParameterSet>(tool_label);
        _signatureToolsVec.push_back(art::make_tool<SignatureToolBase>(tool_pset));
    }

    const auto& clarity_tools_pset = pset.get<fhicl::ParameterSet>("ClarityTools");
    for (const auto& tool_label : clarity_tools_pset.get_pset_names()) {
        auto const& tool_pset = clarity_tools_pset.get<fhicl::ParameterSet>(tool_label);
        _clarityToolsVec.push_back(art::make_tool<ClarityToolBase>(tool_pset));
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

    for (const auto& entry : _pattern) {
        const auto& type = entry.first;
        const auto& sig = entry.second;
        ClarityResult scr(type, sig);

        for (int view = common::TPC_VIEW_U; view != common::N_VIEWS; ++view) {
            for (const auto& clarityTool : _clarityToolsVec) {
                if (!clarityTool->filter(e, sig, type, static_cast<common::PandoraView>(view))) {
                    scr.passes[static_cast<size_t>(view)] = false;
                    break; 
                }
            }
        }
        _clarity_results.push_back(std::move(scr));
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
    if (!_pattern_found) 
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

bool EventClassifier::passesPlaneClarity(const Signature& sig, common::PandoraView view) const {
    auto it = std::find_if(_clarity_results.begin(), _clarity_results.end(),
                           [&sig](const ClarityResult& scr) {
                               return scr.signature == sig;
                           });
    if (it != _clarity_results.end())
        return it->passes[static_cast<size_t>(view)];
    return false;
}

bool EventClassifier::passesThreePlaneClarity(const Signature& sig) const {
    auto it = std::find_if(_clarity_results.begin(), _clarity_results.end(),
                           [&sig](const ClarityResult& scr) {
                               return scr.signature == sig;
                           });
    if (it != _clarity_results.end())
        return std::all_of(it->passes.begin(), it->passes.end(), [](bool passed) { return passed; });
    return false;
}

bool EventClassifier::passClarityFilter() const {
    if (_clarity_results.empty()) 
        return false;
    
    return std::all_of(_clarity_results.begin(), _clarity_results.end(),
        [](const ClarityResult& cr) {
            return std::all_of(cr.passes.begin(), cr.passes.end(),
                [](bool passed) { return passed; });
        });
}

} 

#endif // EVENTCLASSIFIER_H
