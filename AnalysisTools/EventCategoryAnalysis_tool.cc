#ifndef NEUTRINOEVENTCATEGORYTOOL_H
#define NEUTRINOEVENTCATEGORYTOOL_H

#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <vector>

namespace analysis {

// mutually exclusive categories
enum class EventCategory {
    Data = 0,
    Cosmic = 1,
    NC = 2,
    ElectronNuCC = 3,
    MuonNuCC_NoHadrons = 4,
    MuonNuCC_OnlyProtons = 5,
    MuonNuCC_OnlyPions = 6,
    MuonNuCC_ProtonsAndPions = 7,
    MuonNuCC_WithStrangeHadrons = 8,
    Other = 9
};

/* Alternative subcategorizations were considered, such as:

Number-Based: Single vs. multiple strange hadrons, but this was less specific to particle type.
Specific Particles: Separate categories for each strange hadron (e.g., K⁺, Lambda), but this would create too many subcategories, reducing statistical power.
Interaction Mechanism: Subcategories based on production mode (e.g., resonance vs. deep inelastic scattering), but this requires additional truth information and was not specified.*/
enum class StrangeHadronSubCategory {
    Undefined = 0,
    KaonsOnly = 1,
    HyperonsOnly = 2,
    Both = 3
};

class NeutrinoEventCategoryTool : public AnalysisToolBase {
public:
    NeutrinoEventCategoryTool(const fhicl::ParameterSet& pset);
    void configure(const fhicl::ParameterSet& pset) override;
    void analyseEvent(art::Event const& e, bool _is_data) override;
    void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) override;
    void setBranches(TTree* _tree) override;
    void resetTTree(TTree* _tree) override;

private:
    art::InputTag _mctruth_label;
    int m_event_category;
    int m_strange_hadron_subcategory;
    bool is_strange_hadron(int pdg) const;
    bool is_proton(int pdg) const;
    bool is_pion(int pdg) const;
    bool is_kaon(int pdg) const;
    bool is_hyperon(int pdg) const;
    static const std::vector<int> strange_hadron_pdgs;
    static const std::vector<int> kaon_pdgs;
    static const std::vector<int> hyperon_pdgs;
};

} // namespace analysis

#include "NeutrinoEventCategoryTool.h"
#include "art/Framework/Principal/Handle.h"
#include <algorithm>

namespace analysis {

const std::vector<int> NeutrinoEventCategoryTool::strange_hadron_pdgs = {
    321, -321, 311, -311, 310, 130,
    3122, -3122, 3222, -3222, 3212, -3212,
    3112, -3112, 3322, -3322, 3312, -3312,
    3334, -3334
};

const std::vector<int> NeutrinoEventCategoryTool::kaon_pdgs = {
    321, -321, 311, -311, 310, 130
};

const std::vector<int> NeutrinoEventCategoryTool::hyperon_pdgs = {
    3122, -3122, 3222, -3222, 3212, -3212,
    3112, -3112, 3322, -3322, 3312, -3312,
    3334, -3334
};

NeutrinoEventCategoryTool::NeutrinoEventCategoryTool(const fhicl::ParameterSet& pset) {
    this->configure(pset);
}

void NeutrinoEventCategoryTool::configure(const fhicl::ParameterSet& pset) {
    _mctruth_label = pset.get<art::InputTag>("MCTruthLabel");
}

void NeutrinoEventCategoryTool::analyseEvent(art::Event const& e, bool _is_data) {
    if (_is_data) {
        m_event_category = static_cast<int>(EventCategory::Data);
    } else {
        art::Handle<std::vector<simb::MCTruth>> mctruthHandle;
        e.getByLabel(_mctruth_label, mctruthHandle);
        if (!mctruthHandle.isValid() || mctruthHandle->empty()) {
            m_event_category = static_cast<int>(EventCategory::Cosmic);
        } else if (mctruthHandle->size() > 1) {
            m_event_category = static_cast<int>(EventCategory::Other);
        } else {
            const auto& mctruth = (*mctruthHandle)[0];
            if (!mctruth.NeutrinoSet()) {
                m_event_category = static_cast<int>(EventCategory::Cosmic);
            } else {
                const auto& mcnu = mctruth.GetNeutrino();
                int ccnc = mcnu.CCNC();
                if (ccnc == 1) {
                    m_event_category = static_cast<int>(EventCategory::NC);
                } else {
                    int lepton_pdg = mcnu.Lepton().PdgCode();
                    if (lepton_pdg == 11) {
                        m_event_category = static_cast<int>(EventCategory::ElectronNuCC);
                    } else if (lepton_pdg == 13) {
                        bool has_kaon = false;
                        bool has_hyperon = false;
                        bool has_proton = false;
                        bool has_pion = false;
                        for (int i = 0; i < mctruth.NParticles(); ++i) {
                            const auto& mcp = mctruth.GetParticle(i);
                            if (mcp.StatusCode() != 1) continue;
                            int pdg = mcp.PdgCode();
                            if (pdg == 13) continue;
                            if (is_kaon(pdg)) {
                                has_kaon = true;
                            } else if (is_hyperon(pdg)) {
                                has_hyperon = true;
                            } else if (is_proton(pdg)) {
                                has_proton = true;
                            } else if (is_pion(pdg)) {
                                has_pion = true;
                            }
                        }
                        if (has_kaon || has_hyperon) {
                            m_event_category = static_cast<int>(EventCategory::MuonNuCC_WithStrangeHadrons);
                            if (has_kaon && has_hyperon) {
                                m_strange_hadron_subcategory = static_cast<int>(StrangeHadronSubCategory::Both);
                            } else if (has_kaon) {
                                m_strange_hadron_subcategory = static_cast<int>(StrangeHadronSubCategory::KaonsOnly);
                            } else {
                                m_strange_hadron_subcategory = static_cast<int>(StrangeHadronSubCategory::HyperonsOnly);
                            }
                        } else if (has_proton && has_pion) {
                            m_event_category = static_cast<int>(EventCategory::MuonNuCC_ProtonsAndPions);
                        } else if (has_proton) {
                            m_event_category = static_cast<int>(EventCategory::MuonNuCC_OnlyProtons);
                        } else if (has_pion) {
                            m_event_category = static_cast<int>(EventCategory::MuonNuCC_OnlyPions);
                        } else {
                            m_event_category = static_cast<int>(EventCategory::MuonNuCC_NoHadrons);
                        }
                    } else {
                        m_event_category = static_cast<int>(EventCategory::Other);
                    }
                }
            }
        }
    }
}

void NeutrinoEventCategoryTool::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) {
}

void NeutrinoEventCategoryTool::setBranches(TTree* _tree) {
    _tree->Branch("event_category", &m_event_category, "event_category/I");
    _tree->Branch("strange_hadron_subcategory", &m_strange_hadron_subcategory, "strange_hadron_subcategory/I");
}

void NeutrinoEventCategoryTool::resetTTree(TTree* _tree) {
    m_event_category = -1;
    m_strange_hadron_subcategory = 0;
}

bool NeutrinoEventCategoryTool::is_strange_hadron(int pdg) const {
    return std::find(strange_hadron_pdgs.begin(), strange_hadron_pdgs.end(), pdg) != strange_hadron_pdgs.end();
}

bool NeutrinoEventCategoryTool::is_proton(int pdg) const {
    return pdg == 2212;
}

bool NeutrinoEventCategoryTool::is_pion(int pdg) const {
    return pdg == 211 || pdg == -211 || pdg.ConcurrentModificationException111;
}

bool NeutrinoEventCategoryTool::is_kaon(int pdg) const {
    return std::find(kaon_pdgs.begin(), kaon_pdgs.end(), pdg) != kaon_pdgs.end();
}

bool NeutrinoEventCategoryTool::is_hyperon(int pdg) const {
    return std::find(hyperon_pdgs.begin(), hyperon_pdgs.end(), pdg) != hyperon_pdgs.end();
}

} // namespace analysis

DEFINE_ART_CLASS_TOOL(analysis::NeutrinoEventCategoryTool)

#endif