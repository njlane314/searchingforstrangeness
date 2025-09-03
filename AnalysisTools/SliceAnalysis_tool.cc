#ifndef ANALYSIS_SLICEPURCOMPL_CXX
#define ANALYSIS_SLICEPURCOMPL_CXX

#include <iostream>
#include "AnalysisToolBase.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "Common/BacktrackingFuncs.h"

namespace analysis {

class SliceAnalysis : public AnalysisToolBase {
public:
    SliceAnalysis(const fhicl::ParameterSet &pset);
    ~SliceAnalysis() {}
    void configure(fhicl::ParameterSet const &pset);
    void analyseEvent(art::Event const &e, bool is_data) override;
    void analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected) override;
    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

private:
    art::InputTag fCLSproducer;
    art::InputTag fSLCproducer;
    art::InputTag fMCTproducer;
    art::InputTag fMCPproducer;
    art::InputTag fHITproducer;
    art::InputTag fBKTproducer;
    art::InputTag fOrigHITproducer;
    art::InputTag fOrigBKTproducer;

    std::vector<size_t> muon_ids;
    std::vector<size_t> electron_ids;
    std::vector<size_t> proton_ids;
    std::vector<size_t> charged_pion_ids;
    std::vector<size_t> neutral_pion_ids;
    std::vector<size_t> neutron_ids;
    std::vector<size_t> gamma_ids;
    std::vector<size_t> charged_kaon_ids;
    std::vector<size_t> neutral_kaon_ids;
    std::vector<size_t> lambda_ids;
    std::vector<size_t> charged_sigma_ids;
    std::vector<size_t> sigma_zero_ids;
    std::vector<size_t> other_ids;

    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;

    int original_event_neutrino_hits;
    int event_neutrino_hits;
    int event_muon_hits;
    int event_electron_hits;
    int event_proton_hits;
    int event_charged_pion_hits;
    int event_neutral_pion_hits;
    int event_neutron_hits;
    int event_gamma_hits;
    int event_other_hits;
    int event_charged_kaon_hits;
    int event_neutral_kaon_hits;
    int event_lambda_hits;
    int event_charged_sigma_hits;
    int event_sigma_zero_hits;
    int event_cosmic_hits;

    int slice_neutrino_hits;
    int slice_muon_hits;
    int slice_electron_hits;
    int slice_proton_hits;
    int slice_charged_pion_hits;
    int slice_neutral_pion_hits;
    int slice_neutron_hits;
    int slice_gamma_hits;
    int slice_other_hits;
    int slice_charged_kaon_hits;
    int slice_neutral_kaon_hits;
    int slice_lambda_hits;
    int slice_charged_sigma_hits;
    int slice_sigma_zero_hits;
    int slice_cosmic_hits;

    std::vector<int> pfp_neutrino_hits;
    std::vector<int> pfp_muon_hits;
    std::vector<int> pfp_electron_hits;
    std::vector<int> pfp_proton_hits;
    std::vector<int> pfp_charged_pion_hits;
    std::vector<int> pfp_neutral_pion_hits;
    std::vector<int> pfp_neutron_hits;
    std::vector<int> pfp_gamma_hits;
    std::vector<int> pfp_other_hits;
    std::vector<int> pfp_charged_kaon_hits;
    std::vector<int> pfp_neutral_kaon_hits;
    std::vector<int> pfp_lambda_hits;
    std::vector<int> pfp_charged_sigma_hits;
    std::vector<int> pfp_sigma_zero_hits;
    std::vector<int> pfp_cosmic_hits;

    float neutrino_completeness_from_pfp;
    float neutrino_purity_from_pfp;

    enum class PrimaryParticleLabel {
        Empty,
        Cosmic,
        Muon,
        Electron,
        Proton,
        ChargedPion,
        NeutralPion,
        Neutron,
        Gamma,
        ChargedKaon,
        NeutralKaon,
        Lambda,
        ChargedSigma,
        SigmaZero,
        Other
    };

    void incrementCounts(const std::vector<art::Ptr<recob::Hit>> &hits,
                         const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart,
                         const std::vector<size_t> &muon_ids, const std::vector<size_t> &electron_ids,
                         const std::vector<size_t> &proton_ids, const std::vector<size_t> &charged_pion_ids,
                         const std::vector<size_t> &neutral_pion_ids, const std::vector<size_t> &neutron_ids,
                         const std::vector<size_t> &gamma_ids, const std::vector<size_t> &other_ids,
                         const std::vector<size_t> &charged_kaon_ids, const std::vector<size_t> &neutral_kaon_ids,
                         const std::vector<size_t> &lambda_ids, const std::vector<size_t> &charged_sigma_ids,
                         const std::vector<size_t> &sigma_zero_ids,
                         int &muon_hits, int &electron_hits, int &proton_hits, int &charged_pion_hits,
                         int &neutral_pion_hits, int &neutron_hits, int &gamma_hits, int &other_hits,
                         int &charged_kaon_hits, int &neutral_kaon_hits, int &lambda_hits,
                         int &charged_sigma_hits, int &sigma_zero_hits, int &cosmic_hits) const;

    void assignLabelToProgenyRecursively(size_t particle_index,
                                         const std::vector<simb::MCParticle> &particles,
                                         std::vector<PrimaryParticleLabel> &particle_labels,
                                         const std::unordered_map<int, size_t> &track_id_to_index,
                                         PrimaryParticleLabel primary_label_to_assign) const;

    std::vector<PrimaryParticleLabel> classifyParticles(const art::Event &event) const;

    PrimaryParticleLabel getPrimaryLabel(int pdg) const;
};

SliceAnalysis::SliceAnalysis(const fhicl::ParameterSet &p) {
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fSLCproducer = p.get<art::InputTag>("SLCproducer");
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fBKTproducer = p.get<art::InputTag>("BKTproducer");
    fOrigHITproducer = p.get<art::InputTag>("OrigHITproducer");
    fOrigBKTproducer = p.get<art::InputTag>("OrigBKTproducer");
}

void SliceAnalysis::configure(fhicl::ParameterSet const &p) {}

void SliceAnalysis::analyseEvent(art::Event const &e, bool is_data) {
    if (is_data) return;
    art::ValidHandle<std::vector<simb::MCTruth>> inputMCTruth = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    if (inputMCTruth->empty()) return;
    auto particle_labels = this->classifyParticles(e);
    art::ValidHandle<std::vector<simb::MCParticle>> inputMCParticle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
    muon_ids.clear();
    electron_ids.clear();
    proton_ids.clear();
    charged_pion_ids.clear();
    neutral_pion_ids.clear();
    neutron_ids.clear();
    gamma_ids.clear();
    other_ids.clear();
    charged_kaon_ids.clear();
    neutral_kaon_ids.clear();
    lambda_ids.clear();
    charged_sigma_ids.clear();
    sigma_zero_ids.clear();
    for (size_t i = 0; i < inputMCParticle->size(); ++i) {
        const auto &mcp = inputMCParticle->at(i);
        if (mcp.StatusCode() == 1) {
            PrimaryParticleLabel label = particle_labels[i];
            switch (label) {
                case PrimaryParticleLabel::Muon: muon_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::Electron: electron_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::Proton: proton_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::ChargedPion: charged_pion_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::NeutralPion: neutral_pion_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::Neutron: neutron_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::Gamma: gamma_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::Other: other_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::ChargedKaon: charged_kaon_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::NeutralKaon: neutral_kaon_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::Lambda: lambda_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::ChargedSigma: charged_sigma_ids.push_back(mcp.TrackId()); break;
                case PrimaryParticleLabel::SigmaZero: sigma_zero_ids.push_back(mcp.TrackId()); break;
                default: break;
            }
        }
    }
    art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
        new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBKTproducer));
    int muon_hits = 0, electron_hits = 0, proton_hits = 0, charged_pion_hits = 0;
    int neutral_pion_hits = 0, neutron_hits = 0, gamma_hits = 0, other_hits = 0;
    int charged_kaon_hits = 0, neutral_kaon_hits = 0, lambda_hits = 0;
    int charged_sigma_hits = 0, sigma_zero_hits = 0, cosmic_hits = 0;
    std::vector<art::Ptr<recob::Hit>> inHitsPtrV;
    for (unsigned int ih = 0; ih < inputHits->size(); ih++)
        inHitsPtrV.push_back({inputHits, ih});
    this->incrementCounts(inHitsPtrV, assocMCPart,
                          muon_ids, electron_ids, proton_ids, charged_pion_ids,
                          neutral_pion_ids, neutron_ids, gamma_ids, other_ids,
                          charged_kaon_ids, neutral_kaon_ids, lambda_ids,
                          charged_sigma_ids, sigma_zero_ids,
                          muon_hits, electron_hits, proton_hits, charged_pion_hits,
                          neutral_pion_hits, neutron_hits, gamma_hits, other_hits,
                          charged_kaon_hits, neutral_kaon_hits, lambda_hits,
                          charged_sigma_hits, sigma_zero_hits, cosmic_hits);
    event_neutrino_hits = muon_hits + electron_hits + proton_hits + charged_pion_hits +
                          neutral_pion_hits + neutron_hits + gamma_hits + other_hits +
                          charged_kaon_hits + neutral_kaon_hits + lambda_hits +
                          charged_sigma_hits + sigma_zero_hits;
    event_muon_hits = muon_hits;
    event_electron_hits = electron_hits;
    event_proton_hits = proton_hits;
    event_charged_pion_hits = charged_pion_hits;
    event_neutral_pion_hits = neutral_pion_hits;
    event_neutron_hits = neutron_hits;
    event_gamma_hits = gamma_hits;
    event_other_hits = other_hits;
    event_charged_kaon_hits = charged_kaon_hits;
    event_neutral_kaon_hits = neutral_kaon_hits;
    event_lambda_hits = lambda_hits;
    event_charged_sigma_hits = charged_sigma_hits;
    event_sigma_zero_hits = sigma_zero_hits;
    event_cosmic_hits = cosmic_hits;
    if (!fOrigHITproducer.empty()) {
        auto hitsOrig = e.getValidHandle<std::vector<recob::Hit>>(fOrigHITproducer);
        auto assocMCPartOrig = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
            new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitsOrig, e, fOrigBKTproducer));
        int orig_nu_hits = 0;
        for (size_t ih = 0; ih < hitsOrig->size(); ih++) {
            const art::Ptr<recob::Hit> hit_ptr(hitsOrig, ih);
            if (assocMCPartOrig->at(hit_ptr.key()).size()) orig_nu_hits++;
        }
        original_event_neutrino_hits = orig_nu_hits;
    } else {
        original_event_neutrino_hits = event_neutrino_hits;
    }
}

void SliceAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected) {
    if (is_data) return;
    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(
        e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
    art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(
        new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));
    neutrino_purity_from_pfp = 0;
    neutrino_completeness_from_pfp = 0;
    int total_hits = 0;
    for (auto pfp : slice_pfp_v) {
        if (pfp->IsPrimary()) {
            auto slice_pxy_v = pfp.get<recob::Slice>();
            if (slice_pxy_v.size() != 1) {
                std::cout << "WRONG!!! n slices = " << slice_pxy_v.size() << " " << __FILE__ << " " << __LINE__ << std::endl;
                return;
            }
            auto slicehits = assocSliceHit->at(slice_pxy_v[0].key());
            int muon_hits = 0, electron_hits = 0, proton_hits = 0, charged_pion_hits = 0;
            int neutral_pion_hits = 0, neutron_hits = 0, gamma_hits = 0, other_hits = 0;
            int charged_kaon_hits = 0, neutral_kaon_hits = 0, lambda_hits = 0;
            int charged_sigma_hits = 0, sigma_zero_hits = 0, cosmic_hits = 0;
            this->incrementCounts(slicehits, assocMCPart,
                                  muon_ids, electron_ids, proton_ids, charged_pion_ids,
                                  neutral_pion_ids, neutron_ids, gamma_ids, other_ids,
                                  charged_kaon_ids, neutral_kaon_ids, lambda_ids,
                                  charged_sigma_ids, sigma_zero_ids,
                                  muon_hits, electron_hits, proton_hits, charged_pion_hits,
                                  neutral_pion_hits, neutron_hits, gamma_hits, other_hits,
                                  charged_kaon_hits, neutral_kaon_hits, lambda_hits,
                                  charged_sigma_hits, sigma_zero_hits, cosmic_hits);
            slice_neutrino_hits = muon_hits + electron_hits + proton_hits + charged_pion_hits +
                                  neutral_pion_hits + neutron_hits + gamma_hits + other_hits +
                                  charged_kaon_hits + neutral_kaon_hits + lambda_hits +
                                  charged_sigma_hits + sigma_zero_hits;
            slice_muon_hits = muon_hits;
            slice_electron_hits = electron_hits;
            slice_proton_hits = proton_hits;
            slice_charged_pion_hits = charged_pion_hits;
            slice_neutral_pion_hits = neutral_pion_hits;
            slice_neutron_hits = neutron_hits;
            slice_gamma_hits = gamma_hits;
            slice_other_hits = other_hits;
            slice_charged_kaon_hits = charged_kaon_hits;
            slice_neutral_kaon_hits = neutral_kaon_hits;
            slice_lambda_hits = lambda_hits;
            slice_charged_sigma_hits = charged_sigma_hits;
            slice_sigma_zero_hits = sigma_zero_hits;
            slice_cosmic_hits = cosmic_hits;
        }
        std::vector<art::Ptr<recob::Hit>> hit_v;
        auto clus_pxy_v = pfp.get<recob::Cluster>();
        if (clus_pxy_v.size() > 0) {
            for (auto ass_clus : clus_pxy_v) {
                const auto &clus = clus_proxy[ass_clus.key()];
                auto clus_hit_v = clus.get<recob::Hit>();
                for (const auto &hit : clus_hit_v) hit_v.push_back(hit);
            }
        }
        total_hits += hit_v.size();
        int muon_hits = 0, electron_hits = 0, proton_hits = 0, charged_pion_hits = 0;
        int neutral_pion_hits = 0, neutron_hits = 0, gamma_hits = 0, other_hits = 0;
        int charged_kaon_hits = 0, neutral_kaon_hits = 0, lambda_hits = 0;
        int charged_sigma_hits = 0, sigma_zero_hits = 0, cosmic_hits = 0;
        this->incrementCounts(hit_v, assocMCPart,
                              muon_ids, electron_ids, proton_ids, charged_pion_ids,
                              neutral_pion_ids, neutron_ids, gamma_ids, other_ids,
                              charged_kaon_ids, neutral_kaon_ids, lambda_ids,
                              charged_sigma_ids, sigma_zero_ids,
                              muon_hits, electron_hits, proton_hits, charged_pion_hits,
                              neutral_pion_hits, neutron_hits, gamma_hits, other_hits,
                              charged_kaon_hits, neutral_kaon_hits, lambda_hits,
                              charged_sigma_hits, sigma_zero_hits, cosmic_hits);
        int pfp_nu_hits = muon_hits + electron_hits + proton_hits + charged_pion_hits +
                          neutral_pion_hits + neutron_hits + gamma_hits + other_hits +
                          charged_kaon_hits + neutral_kaon_hits + lambda_hits +
                          charged_sigma_hits + sigma_zero_hits;
        pfp_neutrino_hits.push_back(pfp_nu_hits);
        pfp_muon_hits.push_back(muon_hits);
        pfp_electron_hits.push_back(electron_hits);
        pfp_proton_hits.push_back(proton_hits);
        pfp_charged_pion_hits.push_back(charged_pion_hits);
        pfp_neutral_pion_hits.push_back(neutral_pion_hits);
        pfp_neutron_hits.push_back(neutron_hits);
        pfp_gamma_hits.push_back(gamma_hits);
        pfp_other_hits.push_back(other_hits);
        pfp_charged_kaon_hits.push_back(charged_kaon_hits);
        pfp_neutral_kaon_hits.push_back(neutral_kaon_hits);
        pfp_lambda_hits.push_back(lambda_hits);
        pfp_charged_sigma_hits.push_back(charged_sigma_hits);
        pfp_sigma_zero_hits.push_back(sigma_zero_hits);
        pfp_cosmic_hits.push_back(cosmic_hits);
        neutrino_completeness_from_pfp += pfp_nu_hits;
    }
    if (total_hits > 0) neutrino_purity_from_pfp = static_cast<float>(neutrino_completeness_from_pfp) / total_hits;
    if (original_event_neutrino_hits > 0) neutrino_completeness_from_pfp /= original_event_neutrino_hits;
}

void SliceAnalysis::setBranches(TTree *_tree) {
    _tree->Branch("original_event_neutrino_hits", &original_event_neutrino_hits, "original_event_neutrino_hits/I");
    _tree->Branch("event_neutrino_hits", &event_neutrino_hits, "event_neutrino_hits/I");
    _tree->Branch("event_muon_hits", &event_muon_hits, "event_muon_hits/I");
    _tree->Branch("event_electron_hits", &event_electron_hits, "event_electron_hits/I");
    _tree->Branch("event_proton_hits", &event_proton_hits, "event_proton_hits/I");
    _tree->Branch("event_charged_pion_hits", &event_charged_pion_hits, "event_charged_pion_hits/I");
    _tree->Branch("event_neutral_pion_hits", &event_neutral_pion_hits, "event_neutral_pion_hits/I");
    _tree->Branch("event_neutron_hits", &event_neutron_hits, "event_neutron_hits/I");
    _tree->Branch("event_gamma_hits", &event_gamma_hits, "event_gamma_hits/I");
    _tree->Branch("event_other_hits", &event_other_hits, "event_other_hits/I");
    _tree->Branch("event_charged_kaon_hits", &event_charged_kaon_hits, "event_charged_kaon_hits/I");
    _tree->Branch("event_neutral_kaon_hits", &event_neutral_kaon_hits, "event_neutral_kaon_hits/I");
    _tree->Branch("event_lambda_hits", &event_lambda_hits, "event_lambda_hits/I");
    _tree->Branch("event_charged_sigma_hits", &event_charged_sigma_hits, "event_charged_sigma_hits/I");
    _tree->Branch("event_sigma_zero_hits", &event_sigma_zero_hits, "event_sigma_zero_hits/I");
    _tree->Branch("event_cosmic_hits", &event_cosmic_hits, "event_cosmic_hits/I");
    _tree->Branch("slice_neutrino_hits", &slice_neutrino_hits, "slice_neutrino_hits/I");
    _tree->Branch("slice_muon_hits", &slice_muon_hits, "slice_muon_hits/I");
    _tree->Branch("slice_electron_hits", &slice_electron_hits, "slice_electron_hits/I");
    _tree->Branch("slice_proton_hits", &slice_proton_hits, "slice_proton_hits/I");
    _tree->Branch("slice_charged_pion_hits", &slice_charged_pion_hits, "slice_charged_pion_hits/I");
    _tree->Branch("slice_neutral_pion_hits", &slice_neutral_pion_hits, "slice_neutral_pion_hits/I");
    _tree->Branch("slice_neutron_hits", &slice_neutron_hits, "slice_neutron_hits/I");
    _tree->Branch("slice_gamma_hits", &slice_gamma_hits, "slice_gamma_hits/I");
    _tree->Branch("slice_other_hits", &slice_other_hits, "slice_other_hits/I");
    _tree->Branch("slice_charged_kaon_hits", &slice_charged_kaon_hits, "slice_charged_kaon_hits/I");
    _tree->Branch("slice_neutral_kaon_hits", &slice_neutral_kaon_hits, "slice_neutral_kaon_hits/I");
    _tree->Branch("slice_lambda_hits", &slice_lambda_hits, "slice_lambda_hits/I");
    _tree->Branch("slice_charged_sigma_hits", &slice_charged_sigma_hits, "slice_charged_sigma_hits/I");
    _tree->Branch("slice_sigma_zero_hits", &slice_sigma_zero_hits, "slice_sigma_zero_hits/I");
    _tree->Branch("slice_cosmic_hits", &slice_cosmic_hits, "slice_cosmic_hits/I");
    _tree->Branch("pfp_neutrino_hits", &pfp_neutrino_hits);
    _tree->Branch("pfp_muon_hits", &pfp_muon_hits);
    _tree->Branch("pfp_electron_hits", &pfp_electron_hits);
    _tree->Branch("pfp_proton_hits", &pfp_proton_hits);
    _tree->Branch("pfp_charged_pion_hits", &pfp_charged_pion_hits);
    _tree->Branch("pfp_neutral_pion_hits", &pfp_neutral_pion_hits);
    _tree->Branch("pfp_neutron_hits", &pfp_neutron_hits);
    _tree->Branch("pfp_gamma_hits", &pfp_gamma_hits);
    _tree->Branch("pfp_other_hits", &pfp_other_hits);
    _tree->Branch("pfp_charged_kaon_hits", &pfp_charged_kaon_hits);
    _tree->Branch("pfp_neutral_kaon_hits", &pfp_neutral_kaon_hits);
    _tree->Branch("pfp_lambda_hits", &pfp_lambda_hits);
    _tree->Branch("pfp_charged_sigma_hits", &pfp_charged_sigma_hits);
    _tree->Branch("pfp_sigma_zero_hits", &pfp_sigma_zero_hits);
    _tree->Branch("pfp_cosmic_hits", &pfp_cosmic_hits);
    _tree->Branch("neutrino_completeness_from_pfp", &neutrino_completeness_from_pfp, "neutrino_completeness_from_pfp/F");
    _tree->Branch("neutrino_purity_from_pfp", &neutrino_purity_from_pfp, "neutrino_purity_from_pfp/F");
}

void SliceAnalysis::resetTTree(TTree *_tree) {
    original_event_neutrino_hits = std::numeric_limits<int>::min();
    event_neutrino_hits = std::numeric_limits<int>::min();
    event_muon_hits = std::numeric_limits<int>::min();
    event_electron_hits = std::numeric_limits<int>::min();
    event_proton_hits = std::numeric_limits<int>::min();
    event_charged_pion_hits = std::numeric_limits<int>::min();
    event_neutral_pion_hits = std::numeric_limits<int>::min();
    event_neutron_hits = std::numeric_limits<int>::min();
    event_gamma_hits = std::numeric_limits<int>::min();
    event_other_hits = std::numeric_limits<int>::min();
    event_charged_kaon_hits = std::numeric_limits<int>::min();
    event_neutral_kaon_hits = std::numeric_limits<int>::min();
    event_lambda_hits = std::numeric_limits<int>::min();
    event_charged_sigma_hits = std::numeric_limits<int>::min();
    event_sigma_zero_hits = std::numeric_limits<int>::min();
    event_cosmic_hits = std::numeric_limits<int>::min();
    slice_neutrino_hits = std::numeric_limits<int>::min();
    slice_muon_hits = std::numeric_limits<int>::min();
    slice_electron_hits = std::numeric_limits<int>::min();
    slice_proton_hits = std::numeric_limits<int>::min();
    slice_charged_pion_hits = std::numeric_limits<int>::min();
    slice_neutral_pion_hits = std::numeric_limits<int>::min();
    slice_neutron_hits = std::numeric_limits<int>::min();
    slice_gamma_hits = std::numeric_limits<int>::min();
    slice_other_hits = std::numeric_limits<int>::min();
    slice_charged_kaon_hits = std::numeric_limits<int>::min();
    slice_neutral_kaon_hits = std::numeric_limits<int>::min();
    slice_lambda_hits = std::numeric_limits<int>::min();
    slice_charged_sigma_hits = std::numeric_limits<int>::min();
    slice_sigma_zero_hits = std::numeric_limits<int>::min();
    slice_cosmic_hits = std::numeric_limits<int>::min();
    pfp_neutrino_hits.clear();
    pfp_muon_hits.clear();
    pfp_electron_hits.clear();
    pfp_proton_hits.clear();
    pfp_charged_pion_hits.clear();
    pfp_neutral_pion_hits.clear();
    pfp_neutron_hits.clear();
    pfp_gamma_hits.clear();
    pfp_other_hits.clear();
    pfp_charged_kaon_hits.clear();
    pfp_neutral_kaon_hits.clear();
    pfp_lambda_hits.clear();
    pfp_charged_sigma_hits.clear();
    pfp_sigma_zero_hits.clear();
    pfp_cosmic_hits.clear();
    neutrino_completeness_from_pfp = std::numeric_limits<float>::min();
    neutrino_purity_from_pfp = std::numeric_limits<float>::min();
}

void SliceAnalysis::incrementCounts(const std::vector<art::Ptr<recob::Hit>> &hits,
                                    const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart,
                                    const std::vector<size_t> &muon_ids, const std::vector<size_t> &electron_ids,
                                    const std::vector<size_t> &proton_ids, const std::vector<size_t> &charged_pion_ids,
                                    const std::vector<size_t> &neutral_pion_ids, const std::vector<size_t> &neutron_ids,
                                    const std::vector<size_t> &gamma_ids, const std::vector<size_t> &other_ids,
                                    const std::vector<size_t> &charged_kaon_ids, const std::vector<size_t> &neutral_kaon_ids,
                                    const std::vector<size_t> &lambda_ids, const std::vector<size_t> &charged_sigma_ids,
                                    const std::vector<size_t> &sigma_zero_ids,
                                    int &muon_hits, int &electron_hits, int &proton_hits, int &charged_pion_hits,
                                    int &neutral_pion_hits, int &neutron_hits, int &gamma_hits, int &other_hits,
                                    int &charged_kaon_hits, int &neutral_kaon_hits, int &lambda_hits,
                                    int &charged_sigma_hits, int &sigma_zero_hits, int &cosmic_hits) const {
    for (const auto &hit : hits) {
        auto assmcp = assocMCPart->at(hit.key());
        auto assmdt = assocMCPart->data(hit.key());
        if (assmcp.empty()) {
            cosmic_hits++;
            continue;
        }
        for (unsigned int ia = 0; ia < assmcp.size(); ++ia) {
            auto mcp = assmcp[ia];
            auto amd = assmdt[ia];
            if (amd->isMaxIDE != 1) continue;
            size_t track_id = mcp->TrackId();
            if (std::find(muon_ids.begin(), muon_ids.end(), track_id) != muon_ids.end()) muon_hits++;
            else if (std::find(electron_ids.begin(), electron_ids.end(), track_id) != electron_ids.end()) electron_hits++;
            else if (std::find(proton_ids.begin(), proton_ids.end(), track_id) != proton_ids.end()) proton_hits++;
            else if (std::find(charged_pion_ids.begin(), charged_pion_ids.end(), track_id) != charged_pion_ids.end()) charged_pion_hits++;
            else if (std::find(neutral_pion_ids.begin(), neutral_pion_ids.end(), track_id) != neutral_pion_ids.end()) neutral_pion_hits++;
            else if (std::find(neutron_ids.begin(), neutron_ids.end(), track_id) != neutron_ids.end()) neutron_hits++;
            else if (std::find(gamma_ids.begin(), gamma_ids.end(), track_id) != gamma_ids.end()) gamma_hits++;
            else if (std::find(other_ids.begin(), other_ids.end(), track_id) != other_ids.end()) other_hits++;
            else if (std::find(charged_kaon_ids.begin(), charged_kaon_ids.end(), track_id) != charged_kaon_ids.end()) charged_kaon_hits++;
            else if (std::find(neutral_kaon_ids.begin(), neutral_kaon_ids.end(), track_id) != neutral_kaon_ids.end()) neutral_kaon_hits++;
            else if (std::find(lambda_ids.begin(), lambda_ids.end(), track_id) != lambda_ids.end()) lambda_hits++;
            else if (std::find(charged_sigma_ids.begin(), charged_sigma_ids.end(), track_id) != charged_sigma_ids.end()) charged_sigma_hits++;
            else if (std::find(sigma_zero_ids.begin(), sigma_zero_ids.end(), track_id) != sigma_zero_ids.end()) sigma_zero_hits++;
        }
    }
}

void SliceAnalysis::assignLabelToProgenyRecursively(size_t particle_index,
                                                   const std::vector<simb::MCParticle> &particles,
                                                   std::vector<PrimaryParticleLabel> &particle_labels,
                                                   const std::unordered_map<int, size_t> &track_id_to_index,
                                                   PrimaryParticleLabel primary_label_to_assign) const {
    if (particle_index >= particles.size() || particle_index >= particle_labels.size()) return;
    particle_labels[particle_index] = primary_label_to_assign;
    const auto &particle = particles[particle_index];
    for (int daughter_idx = 0; daughter_idx < particle.NumberDaughters(); ++daughter_idx) {
        int daughter_track_id = particle.Daughter(daughter_idx);
        auto it = track_id_to_index.find(daughter_track_id);
        if (it != track_id_to_index.end() && it->second < particles.size()) {
            this->assignLabelToProgenyRecursively(it->second, particles, particle_labels, track_id_to_index, primary_label_to_assign);
        }
    }
}

std::vector<SliceAnalysis::PrimaryParticleLabel> SliceAnalysis::classifyParticles(const art::Event &event) const {
    const auto particle_collection_handle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
    const auto &particles = *particle_collection_handle;
    std::unordered_map<int, size_t> track_id_to_vector_index;
    for (size_t i = 0; i < particles.size(); ++i) {
        track_id_to_vector_index[particles[i].TrackId()] = i;
    }
    std::vector<PrimaryParticleLabel> classified_particle_labels(particles.size(), PrimaryParticleLabel::Empty);
    for (size_t i = 0; i < particles.size(); ++i) {
        if (particles[i].Mother() == 0) {
            auto it = track_id_to_vector_index.find(particles[i].TrackId());
            if (it != track_id_to_vector_index.end()) {
                PrimaryParticleLabel initial_label = this->getPrimaryLabel(particles[i].PdgCode());
                this->assignLabelToProgenyRecursively(it->second, particles, classified_particle_labels, track_id_to_vector_index, initial_label);
            }
        }
    }
    return classified_particle_labels;
}

SliceAnalysis::PrimaryParticleLabel SliceAnalysis::getPrimaryLabel(int pdg) const {
    int abs_pdg = std::abs(pdg);
    if (abs_pdg == 13) return PrimaryParticleLabel::Muon;
    else if (abs_pdg == 11) return PrimaryParticleLabel::Electron;
    else if (abs_pdg == 2212) return PrimaryParticleLabel::Proton;
    else if (abs_pdg == 211) return PrimaryParticleLabel::ChargedPion;
    else if (abs_pdg == 111) return PrimaryParticleLabel::NeutralPion;
    else if (abs_pdg == 2112) return PrimaryParticleLabel::Neutron;
    else if (abs_pdg == 22) return PrimaryParticleLabel::Gamma;
    else if (abs_pdg == 321) return PrimaryParticleLabel::ChargedKaon;
    else if (abs_pdg == 311 || abs_pdg == 310 || abs_pdg == 130) return PrimaryParticleLabel::NeutralKaon;
    else if (abs_pdg == 3122) return PrimaryParticleLabel::Lambda;
    else if (abs_pdg == 3222 || abs_pdg == 3112) return PrimaryParticleLabel::ChargedSigma;
    else if (abs_pdg == 3212) return PrimaryParticleLabel::SigmaZero;
    else if (abs_pdg == 3322 || abs_pdg == 3312 || abs_pdg == 3334) return PrimaryParticleLabel::Other;
    else return PrimaryParticleLabel::Other;
}

DEFINE_ART_CLASS_TOOL(SliceAnalysis)

}

#endif
