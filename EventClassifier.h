#ifndef EVENTCLASSIFIER_H
#define EVENTCLASSIFIER_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "ClarityTools/ClarityToolBase.h"
#include "SignatureTools/SignatureToolBase.h"
#include <vector>
#include <string>
#include <limits>
#include <array>
#include <algorithm>
#include <memory>
#include <iostream>

namespace signature {

enum EventType {
    kSignal = 0,
    kBeamNeutrino,
    kCosmicRay,
    kOther
};

inline bool compareSignatures(const Signature& one_sig, const Signature& two_sig) {
    if (one_sig.size() != two_sig.size())
        return false;
    for (size_t i = 0; i < one_sig.size(); ++i) {
        if (one_sig[i]->TrackId() != two_sig[i]->TrackId())
            return false;
    }
    return true;
}

class ClarityEvaluator {
public:
    class Result {
    public:
        Result(SignatureType t, const Signature& s)
            : type(t), signature(s) { 
                pass_flags.fill(true); 
            }
        bool passes(common::PandoraView view) const { 
            return pass_flags[static_cast<size_t>(view)]; 
        }
        bool allPass() const { 
            return std::all_of(pass_flags.begin(), pass_flags.end(), [](bool p){ return p; }); 
        }
        Signature signature;
        SignatureType type;
        std::array<bool, common::N_VIEWS> pass_flags;
    };

    explicit ClarityEvaluator(const std::vector<std::unique_ptr<ClarityToolBase>>& clarityTools)
      : _clarityToolsVec(clarityTools) {}

    Result evaluate(const art::Event& e, SignatureType type, const Signature& sig) const {
        Result result(type, sig);
        for (int view = common::TPC_VIEW_U; view < common::N_VIEWS; ++view) {
            bool view_pass = true;
            for (const auto& tool : _clarityToolsVec) {
                if (!tool->filter(e, sig, type, static_cast<common::PandoraView>(view))) {
                    view_pass = false;
                    break; 
                }
            }
            result.pass_flags[static_cast<size_t>(view)] = view_pass;
        }
        return result;
    }
private:
    const std::vector<std::unique_ptr<ClarityToolBase>>& _clarityToolsVec;
};

class EventClassifier {
public:
    explicit EventClassifier(fhicl::ParameterSet const& pset);
    EventType classifyEvent(const art::Event& e) const;
    const Pattern& getPattern(const art::Event& e) const;
    bool passClarity(const Signature& sig, common::PandoraView view) const;
    bool passClarity() const;
private:
    simb::Origin_t getTruthOrigin(const art::Event& e) const;
    void createPattern(const art::Event& e) const;
    bool isContained(const art::Event& e) const;
    bool isSignal(const art::Event& e) const;
    
    art::InputTag _MCTproducer;
    std::vector<float> _fiducial_offsets;
    std::vector<std::unique_ptr<SignatureToolBase>> _signatureToolsVec;
    std::vector<std::unique_ptr<ClarityToolBase>> _clarityToolsVec;
    mutable Pattern _pattern;
    mutable std::vector<bool> _signature_detectable;
    std::unique_ptr<ClarityEvaluator> _clarityEvaluator;
    mutable std::vector<ClarityEvaluator::Result> _clarity_results;
};

EventClassifier::EventClassifier(fhicl::ParameterSet const& pset)
    : _fiducial_offsets{pset.get<std::vector<float>>("FiducialOffsets", {10,10,10,10,10,10})},
      _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")} {
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
    _clarityEvaluator = std::make_unique<ClarityEvaluator>(_clarityToolsVec);
}

EventType EventClassifier::classifyEvent(const art::Event& e) const {
    createPattern(e);
    if (isSignal(e)) 
        return kSignal;
    simb::Origin_t origin = getTruthOrigin(e);
    switch (origin) {
        case simb::kBeamNeutrino: return kBeamNeutrino;
        case simb::kCosmicRay:    return kCosmicRay;
        default:                  return kOther;
    }
}

void EventClassifier::createPattern(const art::Event& e) const {
    _pattern.clear();
    _signature_detectable.clear();
    _clarity_results.clear();
    for (const auto& tool : _signatureToolsVec) {
        Signature signature;
        bool found = tool->constructSignature(e, signature);
        _pattern.emplace_back(tool->getSignatureType(), signature);
        _signature_detectable.push_back(found ? tool->isDetectable(e, signature) : false);
    }
    for (const auto& entry : _pattern) {
        const auto& type = entry.first;
        const auto& sig  = entry.second;
        _clarity_results.push_back(_clarityEvaluator->evaluate(e, type, sig));
    }
}

const Pattern& EventClassifier::getPattern(const art::Event& e) const {
    createPattern(e);
    return _pattern;
}

simb::Origin_t EventClassifier::getTruthOrigin(const art::Event& e) const {
    art::Handle<std::vector<simb::MCTruth>> truthHandle;
    e.getByLabel(_MCTproducer, truthHandle);
    if (!truthHandle.isValid() || truthHandle->empty())
        return simb::kUnknown;
    return truthHandle->front().Origin();
}

bool EventClassifier::isSignal(const art::Event& e) const {
    createPattern(e);
    bool completePattern = std::all_of(_pattern.begin(), _pattern.end(),
        [](const std::pair<SignatureType, Signature>& entry) {
            return !entry.second.empty();
        });
    return isContained(e) && completePattern &&
           std::all_of(_signature_detectable.begin(), _signature_detectable.end(), [](bool d){ return d; });
}

bool EventClassifier::isContained(const art::Event& e) const {
    art::Handle<std::vector<simb::MCTruth>> truthHandle;
    e.getByLabel(_MCTproducer, truthHandle);
    if (!truthHandle.isValid() || truthHandle->empty())
        return false;
    const simb::MCTruth& truth = truthHandle->front();
    const auto& neutrino = truth.GetNeutrino();
    double point[3] = { neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz() };
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

bool EventClassifier::passClarity(const Signature& sig, common::PandoraView view) const {
    auto it = std::find_if(_clarity_results.begin(), _clarity_results.end(),
        [&sig](const ClarityEvaluator::Result& result) {
            return compareSignatures(result.signature, sig);
        });
    return (it != _clarity_results.end()) ? it->passes(view) : false;
}

bool EventClassifier::passClarity() const {
    if (_clarity_results.empty()) 
        return false;
    return std::all_of(_clarity_results.begin(), _clarity_results.end(),
        [](const ClarityEvaluator::Result& result) {
            return result.allPass();
        });
}

} 

#endif // EVENTCLASSIFIER_H