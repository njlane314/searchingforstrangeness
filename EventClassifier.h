#ifndef EVENTCLASSIFIER_H
#define EVENTCLASSIFIER_H

#include <vector>
#include <memory>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "ClarityTools/ClarityToolBase.h"
#include "SignatureTools/SignatureToolBase.h"

namespace signature 
{
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

    class EventClassifier {
    public:
        explicit EventClassifier(const fhicl::ParameterSet& pset)
            : _fiducialOffsets{pset.get<std::vector<float>>("FiducialOffsets", {10,10,10,10,10,10})},
            _mctProducer{pset.get<art::InputTag>("MCTproducer", "generator")},
            _clarityResultsPtr(nullptr) {
            const auto& sigToolsPset = pset.get<fhicl::ParameterSet>("SignatureTools");
            for (const auto& toolLabel : sigToolsPset.get_pset_names()) {
                auto const& toolPset = sigToolsPset.get<fhicl::ParameterSet>(toolLabel);
                _signatureToolsVec.push_back(art::make_tool<SignatureToolBase>(toolPset));
            }
            const auto& clarityToolsPset = pset.get<fhicl::ParameterSet>("ClarityTools");
            for (const auto& toolLabel : clarityToolsPset.get_pset_names()) {
                auto const& toolPset = clarityToolsPset.get<fhicl::ParameterSet>(toolLabel);
                _clarityToolsVec.push_back(art::make_tool<ClarityToolBase>(toolPset));
            }
            _clarityEvaluator = std::make_unique<ClarityEvaluator>(_clarityToolsVec);
        }

        EventType classifyEvent(const art::Event& e) const {
            createPattern(e);
            if (isSignal(e)) {
                return kSignal;
            }
            simb::Origin_t origin = getTruthOrigin(e);
            switch (origin) {
                case simb::kBeamNeutrino: return kBeamNeutrino;
                case simb::kCosmicRay:    return kCosmicRay;
                default:                  return kOther;
            }
        }

        const Pattern& getPattern(const art::Event& e) const {
            createPattern(e);
            return _pattern;
        }

        bool isSignal(const art::Event& e) const {
            createPattern(e);
            bool completePattern = std::all_of(_pattern.begin(), _pattern.end(),
                [](const std::pair<SignatureType, Signature>& entry) {
                    return !entry.second.empty();
                });
            return isContained(e) && completePattern &&
                std::all_of(_signatureDetectable.begin(), _signatureDetectable.end(), [](bool d){ return d; });
        }

        bool passClarity(const Signature& sig, common::PandoraView view) const {
            if (!_clarityResultsPtr) {
                return false;
            }
            auto it = std::find_if(_clarityResultsPtr->begin(), _clarityResultsPtr->end(),
                [&sig](const ClarityEvaluator::Result& result) {
                    return compareSignatures(result.signature, sig);
                });
            return (it != _clarityResultsPtr->end()) ? it->passes(view) : false;
        }

        bool passClarity() const {
            if (!_clarityResultsPtr || _clarityResultsPtr->empty()) {
                return false;
            }
            return std::all_of(_clarityResultsPtr->begin(), _clarityResultsPtr->end(),
                [](const ClarityEvaluator::Result& result) {
                    return result.allPass();
                });
        }

        const std::vector<ClarityEvaluator::Result>& getClarityResults() const {
            if (!_clarityResultsPtr) {
                static std::vector<ClarityEvaluator::Result> empty;
                return empty;
            }
            return *_clarityResultsPtr;
        }

    private:
        simb::Origin_t getTruthOrigin(const art::Event& e) const {
            art::Handle<std::vector<simb::MCTruth>> truthHandle;
            e.getByLabel(_mctProducer, truthHandle);
            if (!truthHandle.isValid() || truthHandle->empty()) {
                return simb::kUnknown;
            }
            return truthHandle->front().Origin();
        }

        void createPattern(const art::Event& e) const {
            _pattern.clear();
            _signatureDetectable.clear();
            for (const auto& tool : _signatureToolsVec) {
                Signature signature;
                bool found = tool->constructSignature(e, signature);
                _pattern.emplace_back(tool->getSignatureType(), signature);
                _signatureDetectable.push_back(found ? tool->isDetectable(e, signature) : false);
            }
            for (const auto& entry : _pattern) {
                const auto& type = entry.first;
                const auto& sig  = entry.second;
                _clarityEvaluator->evaluate(e, type, sig);
            }
            _clarityResultsPtr = &_clarityEvaluator->getResults();
        }

        bool isContained(const art::Event& e) const {
            art::Handle<std::vector<simb::MCTruth>> truthHandle;
            e.getByLabel(_mctProducer, truthHandle);
            if (!truthHandle.isValid() || truthHandle->empty()) {
                return false;
            }
            const simb::MCTruth& truth = truthHandle->front();
            const auto& neutrino = truth.GetNeutrino();
            double point[3] = {neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz()};
            art::ServiceHandle<geo::Geometry> geo;
            geo::TPCGeo const &tpc = geo->TPC();
            geo::BoxBoundedGeo boundBox = tpc.ActiveBoundingBox();
            std::vector<double> bounds = {boundBox.MinX(), boundBox.MaxX(),
                                        boundBox.MinY(), boundBox.MaxY(),
                                        boundBox.MinZ(), boundBox.MaxZ()};
            bool isX = point[0] > (bounds[0] + _fiducialOffsets[0]) && point[0] < (bounds[1] - _fiducialOffsets[3]);
            bool isY = point[1] > (bounds[2] + _fiducialOffsets[1]) && point[1] < (bounds[3] - _fiducialOffsets[4]);
            bool isZ = point[2] > (bounds[4] + _fiducialOffsets[2]) && point[2] < (bounds[5] - _fiducialOffsets[5]);
            return isX && isY && isZ;
        }

        std::vector<float> _fiducialOffsets;
        art::InputTag _mctProducer;
        std::vector<std::unique_ptr<SignatureToolBase>> _signatureToolsVec;
        std::vector<std::unique_ptr<ClarityToolBase>> _clarityToolsVec;
        mutable Pattern _pattern;
        mutable std::vector<bool> _signatureDetectable;
        std::unique_ptr<ClarityEvaluator> _clarityEvaluator;
        mutable const std::vector<ClarityEvaluator::Result>* _clarityResultsPtr;
    };
}

#endif