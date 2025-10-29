#ifndef ANALYSIS_SLICEPURCOMPL_CXX
#define ANALYSIS_SLICEPURCOMPL_CXX

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "AnalysisToolBase.h"
#include "Common/BacktrackingUtilities.h"

#include <iostream>
#include <limits>
#include <unordered_map>

namespace analysis {

class SliceAnalysis : public AnalysisToolBase {
public:
    SliceAnalysis(const fhicl::ParameterSet &pset);
    ~SliceAnalysis() {}
    void configure(fhicl::ParameterSet const &pset);
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool is_data, bool is_selected) override;
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

    struct HitCategoryCounts {
        int muon = 0;
        int electron = 0;
        int proton = 0;
        int charged_pion = 0;
        int neutral_pion = 0;
        int neutron = 0;
        int gamma = 0;
        int other = 0;
        int charged_kaon = 0;
        int neutral_kaon = 0;
        int lambda_baryon = 0;
        int charged_sigma = 0;
        int sigma_zero = 0;
        int cosmic = 0;
    };

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
                         HitCategoryCounts &out) const;
    static int neutrinoTotal(const HitCategoryCounts &counts);

    void propagateLabel(size_t particle_index,
                        const std::vector<simb::MCParticle> &particles,
                        std::vector<PrimaryParticleLabel> &particle_labels,
                        const std::unordered_map<int, size_t> &track_id_to_index,
                        PrimaryParticleLabel primary_label_to_assign) const;

    std::vector<PrimaryParticleLabel> classifyParticles(const art::Event &event) const;

    PrimaryParticleLabel getPrimaryLabel(int pdg) const;

    std::unordered_map<size_t, PrimaryParticleLabel> fTrackLabelById;
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

void SliceAnalysis::analyseEvent(const art::Event &event, bool is_data) {
    if (is_data) return;
    art::ValidHandle<std::vector<simb::MCTruth>> inputMCTruth = event.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    if (inputMCTruth->empty()) return;
    auto particle_labels = this->classifyParticles(event);
    art::ValidHandle<std::vector<simb::MCParticle>> inputMCParticle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
    muon_ids.clear();
    electron_ids.clear();
    fTrackLabelById.clear();
    fTrackLabelById.reserve(inputMCParticle->size());
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
        if (mcp.StatusCode() != 1) continue;
        const PrimaryParticleLabel label = particle_labels[i];
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
        fTrackLabelById.emplace(static_cast<size_t>(mcp.TrackId()), label);
    }
    art::ValidHandle<std::vector<recob::Hit>> inputHits = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
        new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, event, fBKTproducer));
    std::vector<art::Ptr<recob::Hit>> inHitsPtrV;
    for (unsigned int ih = 0; ih < inputHits->size(); ih++)
        inHitsPtrV.push_back({inputHits, ih});
    HitCategoryCounts event_counts;
    this->incrementCounts(inHitsPtrV, event_counts);
    event_neutrino_hits = neutrinoTotal(event_counts);
    event_muon_hits = event_counts.muon;
    event_electron_hits = event_counts.electron;
    event_proton_hits = event_counts.proton;
    event_charged_pion_hits = event_counts.charged_pion;
    event_neutral_pion_hits = event_counts.neutral_pion;
    event_neutron_hits = event_counts.neutron;
    event_gamma_hits = event_counts.gamma;
    event_other_hits = event_counts.other;
    event_charged_kaon_hits = event_counts.charged_kaon;
    event_neutral_kaon_hits = event_counts.neutral_kaon;
    event_lambda_hits = event_counts.lambda_baryon;
    event_charged_sigma_hits = event_counts.charged_sigma;
    event_sigma_zero_hits = event_counts.sigma_zero;
    event_cosmic_hits = event_counts.cosmic;
    if (!fOrigHITproducer.empty()) {
        auto hitsOrig = event.getValidHandle<std::vector<recob::Hit>>(fOrigHITproducer);
        auto assocMCPartOrig = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
            new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitsOrig, event, fOrigBKTproducer));
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

void SliceAnalysis::analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool is_data, bool is_selected) {
    if (is_data) return;
    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(
        event, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
    art::ValidHandle<std::vector<recob::Slice>> inputSlice = event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(
        new art::FindManyP<recob::Hit>(inputSlice, event, fSLCproducer));
    neutrino_purity_from_pfp = 0;
    neutrino_completeness_from_pfp = 0;
    int total_hits = 0;
    for (auto pfp : slice_pfp_vec) {
        if (pfp->IsPrimary()) {
            auto slice_pxy_v = pfp.get<recob::Slice>();
            if (slice_pxy_v.size() != 1) {
                std::cout << "WRONG!!! n slices = " << slice_pxy_v.size() << " " << __FILE__ << " " << __LINE__ << std::endl;
                return;
            }
            auto slicehits = assocSliceHit->at(slice_pxy_v[0].key());
            HitCategoryCounts slice_counts;
            this->incrementCounts(slicehits, slice_counts);
            slice_neutrino_hits = neutrinoTotal(slice_counts);
            slice_muon_hits = slice_counts.muon;
            slice_electron_hits = slice_counts.electron;
            slice_proton_hits = slice_counts.proton;
            slice_charged_pion_hits = slice_counts.charged_pion;
            slice_neutral_pion_hits = slice_counts.neutral_pion;
            slice_neutron_hits = slice_counts.neutron;
            slice_gamma_hits = slice_counts.gamma;
            slice_other_hits = slice_counts.other;
            slice_charged_kaon_hits = slice_counts.charged_kaon;
            slice_neutral_kaon_hits = slice_counts.neutral_kaon;
            slice_lambda_hits = slice_counts.lambda_baryon;
            slice_charged_sigma_hits = slice_counts.charged_sigma;
            slice_sigma_zero_hits = slice_counts.sigma_zero;
            slice_cosmic_hits = slice_counts.cosmic;
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
        HitCategoryCounts pfp_counts;
        this->incrementCounts(hit_v, pfp_counts);
        int pfp_nu_hits = neutrinoTotal(pfp_counts);
        pfp_neutrino_hits.push_back(pfp_nu_hits);
        pfp_muon_hits.push_back(pfp_counts.muon);
        pfp_electron_hits.push_back(pfp_counts.electron);
        pfp_proton_hits.push_back(pfp_counts.proton);
        pfp_charged_pion_hits.push_back(pfp_counts.charged_pion);
        pfp_neutral_pion_hits.push_back(pfp_counts.neutral_pion);
        pfp_neutron_hits.push_back(pfp_counts.neutron);
        pfp_gamma_hits.push_back(pfp_counts.gamma);
        pfp_other_hits.push_back(pfp_counts.other);
        pfp_charged_kaon_hits.push_back(pfp_counts.charged_kaon);
        pfp_neutral_kaon_hits.push_back(pfp_counts.neutral_kaon);
        pfp_lambda_hits.push_back(pfp_counts.lambda_baryon);
        pfp_charged_sigma_hits.push_back(pfp_counts.charged_sigma);
        pfp_sigma_zero_hits.push_back(pfp_counts.sigma_zero);
        pfp_cosmic_hits.push_back(pfp_counts.cosmic);
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
    original_event_neutrino_hits = -1;
    event_neutrino_hits = -1;
    event_muon_hits = -1;
    event_electron_hits = -1;
    event_proton_hits = -1;
    event_charged_pion_hits = -1;
    event_neutral_pion_hits = -1;
    event_neutron_hits = -1;
    event_gamma_hits = -1;
    event_other_hits = -1;
    event_charged_kaon_hits = -1;
    event_neutral_kaon_hits = -1;
    event_lambda_hits = -1;
    event_charged_sigma_hits = -1;
    event_sigma_zero_hits = -1;
    event_cosmic_hits = -1;
    slice_neutrino_hits = -1;
    slice_muon_hits = -1;
    slice_electron_hits = -1;
    slice_proton_hits = -1;
    slice_charged_pion_hits = -1;
    slice_neutral_pion_hits = -1;
    slice_neutron_hits = -1;
    slice_gamma_hits = -1;
    slice_other_hits = -1;
    slice_charged_kaon_hits = -1;
    slice_neutral_kaon_hits = -1;
    slice_lambda_hits = -1;
    slice_charged_sigma_hits = -1;
    slice_sigma_zero_hits = -1;
    slice_cosmic_hits = -1;
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
    neutrino_completeness_from_pfp = std::numeric_limits<float>::quiet_NaN();
    neutrino_purity_from_pfp = std::numeric_limits<float>::quiet_NaN();
}

void SliceAnalysis::incrementCounts(const std::vector<art::Ptr<recob::Hit>> &hits,
                                    HitCategoryCounts &out) const {
    if (!assocMCPart) return;
    for (const auto &hit : hits) {
        const auto &assmcp = assocMCPart->at(hit.key());
        const auto &assmdt = assocMCPart->data(hit.key());
        if (assmcp.empty()) {
            ++out.cosmic;
            continue;
        }
        for (size_t ia = 0; ia < assmcp.size(); ++ia) {
            auto mcp = assmcp[ia];
            auto amd = assmdt[ia];
            if (amd->isMaxIDE != 1) continue;
            const size_t track_id = static_cast<size_t>(mcp->TrackId());
            auto it = fTrackLabelById.find(track_id);
            if (it == fTrackLabelById.end()) continue;
            switch (it->second) {
                case PrimaryParticleLabel::Muon: ++out.muon; break;
                case PrimaryParticleLabel::Electron: ++out.electron; break;
                case PrimaryParticleLabel::Proton: ++out.proton; break;
                case PrimaryParticleLabel::ChargedPion: ++out.charged_pion; break;
                case PrimaryParticleLabel::NeutralPion: ++out.neutral_pion; break;
                case PrimaryParticleLabel::Neutron: ++out.neutron; break;
                case PrimaryParticleLabel::Gamma: ++out.gamma; break;
                case PrimaryParticleLabel::ChargedKaon: ++out.charged_kaon; break;
                case PrimaryParticleLabel::NeutralKaon: ++out.neutral_kaon; break;
                case PrimaryParticleLabel::Lambda: ++out.lambda_baryon; break;
                case PrimaryParticleLabel::ChargedSigma: ++out.charged_sigma; break;
                case PrimaryParticleLabel::SigmaZero: ++out.sigma_zero; break;
                case PrimaryParticleLabel::Other: ++out.other; break;
                default: break;
            }
        }
    }
}

int SliceAnalysis::neutrinoTotal(const HitCategoryCounts &counts) {
    return counts.muon + counts.electron + counts.proton + counts.charged_pion +
           counts.neutral_pion + counts.neutron + counts.gamma + counts.other +
           counts.charged_kaon + counts.neutral_kaon + counts.lambda_baryon +
           counts.charged_sigma + counts.sigma_zero;
}

void SliceAnalysis::propagateLabel(size_t particle_index,
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
            this->propagateLabel(it->second, particles, particle_labels, track_id_to_index, primary_label_to_assign);
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
                this->propagateLabel(it->second, particles, classified_particle_labels, track_id_to_vector_index, initial_label);
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
