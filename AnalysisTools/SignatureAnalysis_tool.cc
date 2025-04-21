#ifndef EVENTSIGNATUREANALYSISTOOL_H
#define EVENTSIGNATUREANALYSISTOOL_H

#include "AnalysisToolBase.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "../SignatureTools/SignatureToolBase.h"
#include "../ClarityTools/ClarityToolBase.h"
#include "../CommonFunctions/Types.h"
#include <vector>
#include <memory>

namespace analysis {

class EventSignatureAnalysisTool : public AnalysisToolBase {
public:
    EventSignatureAnalysisTool(const fhicl::ParameterSet& pset);
    void configure(const fhicl::ParameterSet& pset) override;
    void analyseEvent(art::Event const& e, bool _is_data) override;
    void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) override;
    void setBranches(TTree* _tree) override;
    void resetTTree(TTree* _tree) override;

private:
    art::InputTag _mctruth_label;
    std::unique_ptr<signature::SignatureToolBase> _signature_tool;
    std::unique_ptr<signature::ClarityToolBase> _clarity_tool;
    int m_signature_type;
    double m_exclusivity_ratio;
    double m_hit_exclusivity_fraction;
    bool m_passes_clarity;
};

} // namespace analysis

#endif // EVENTSIGNATUREANALYSISTOOL_H


#include "EventSignatureAnalysisTool.h"
#include "art/Framework/Principal/Handle.h"
#include <iostream>

namespace analysis {

EventSignatureAnalysisTool::EventSignatureAnalysisTool(const fhicl::ParameterSet& pset) {
    this->configure(pset);
}

void EventSignatureAnalysisTool::configure(const fhicl::ParameterSet& pset) {
    _mctruth_label = pset.get<art::InputTag>("MCTruthLabel");
    auto signature_tool_pset = pset.get<fhicl::ParameterSet>("SignatureTool");
    _signature_tool = art::make_tool<signature::SignatureToolBase>(signature_tool_pset);
    auto clarity_tool_pset = pset.get<fhicl::ParameterSet>("ClarityTool");
    _clarity_tool = art::make_tool<signature::ClarityToolBase>(clarity_tool_pset);
}

void EventSignatureAnalysisTool::analyseEvent(art::Event const& e, bool _is_data) {
    if (_is_data) {
        m_signature_type = -1;
        m_exclusivity_ratio = -1.0;
        m_hit_exclusivity_fraction = -1.0;
        m_passes_clarity = false;
        return;
    }

    signature::Signature sig;
    bool has_signature = _signature_tool->constructSignature(e, sig);
    m_signature_type = static_cast<int>(_signature_tool->getSignatureType());
    if (!has_signature) {
        m_exclusivity_ratio = -1.0;
        m_hit_exclusivity_fraction = -1.0;
        m_passes_clarity = false;
        return;
    }

    common::PandoraView view = common::TPC_VIEW_U;
    m_passes_clarity = _clarity_tool->filter(e, sig, _signature_tool->getSignatureType(), view);
    auto metrics = _clarity_tool->getMetrics();
    auto exclusivity_metrics = dynamic_cast<signature::ExclusivityMetrics*>(metrics.get());
    if (exclusivity_metrics) {
        m_exclusivity_ratio = exclusivity_metrics->exclusivity_ratio;
        m_hit_exclusivity_fraction = exclusivity_metrics->hit_exclusivity_fraction;
    } else {
        m_exclusivity_ratio = -1.0;
        m_hit_exclusivity_fraction = -1.0;
    }
}

void EventSignatureAnalysisTool::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) {
    // No slice-level analysis needed for this tool
}

void EventSignatureAnalysisTool::setBranches(TTree* _tree) {
    _tree->Branch("signature_type", &m_signature_type, "signature_type/I");
    _tree->Branch("exclusivity_ratio", &m_exclusivity_ratio, "exclusivity_ratio/D");
    _tree->Branch("hit_exclusivity_fraction", &m_hit_exclusivity_fraction, "hit_exclusivity_fraction/D");
    _tree->Branch("passes_clarity", &m_passes_clarity, "passes_clarity/O");
}

void EventSignatureAnalysisTool::resetTTree(TTree* _tree) {
    m_signature_type = -1;
    m_exclusivity_ratio = -1.0;
    m_hit_exclusivity_fraction = -1.0;
    m_passes_clarity = false;
}

} // namespace analysis

DEFINE_ART_CLASS_TOOL(analysis::EventSignatureAnalysisTool)