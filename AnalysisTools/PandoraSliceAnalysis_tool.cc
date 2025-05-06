#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "../CommonDefs/Types.h"
#include "TTree.h"
#include <vector>
#include "art/Framework/Principal/Handle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "../CommonDefs/BacktrackingFuncs.h"

namespace analysis 
{
    class PandoraSliceAnalysis : public AnalysisToolBase {
    public:
        explicit PandoraSliceAnalysis(fhicl::ParameterSet const& p);
        virtual ~PandoraSliceAnalysis() = default;

        void configure(const fhicl::ParameterSet& p) override;

        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) override;

        void analyseEvent(art::Event const& e, bool _is_data) override;
        
        void setBranches(TTree* _tree) override;

        void resetTTree(TTree* _tree) override;

    private:
        art::InputTag fPFPproducer;
        art::InputTag fCLSproducer;
        art::InputTag fHproducer;
        art::InputTag fBacktrackTag;
        art::InputTag fMCRproducer;
        
        float slice_nu_score;
        float slice_purity;
        float slice_completeness;
    };

    PandoraSliceAnalysis::PandoraSliceAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void PandoraSliceAnalysis::configure(const fhicl::ParameterSet& p) {
        fPFPproducer = p.get<art::InputTag>("PFPproducer", "pandora");
        fCLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
        fHproducer = p.get<art::InputTag>("Hproducer", "gaushit");
        fBacktrackTag = p.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
        fMCRproducer = p.get<art::InputTag>("MCRproducer", "largeant");
    }

    void PandoraSliceAnalysis::analyseEvent(art::Event const& e, bool _is_data) {}

    void PandoraSliceAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) {
        slice_nu_score = -1.0;
        slice_purity = -1.0;
        slice_completeness = -1.0;

        for (const auto& pfp : slice_pfp_v) {
            if (pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 12 || std::abs(pfp->PdgCode()) == 14)) {
                auto metadata_pxy_v = pfp.get<larpandoraobj::PFParticleMetadata>();
                if (metadata_pxy_v.size() > 0) {
                    const art::Ptr<larpandoraobj::PFParticleMetadata>& pfParticleMetadata = metadata_pxy_v.at(0);
                    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
                    if (pfParticlePropertiesMap.find("NuScore") != pfParticlePropertiesMap.end()) {
                        slice_nu_score = pfParticlePropertiesMap["NuScore"];
                    }
                }
                break;
            }
        }

        if (!_is_data) {
            common::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                                         proxy::withAssociated<recob::Hit>(fCLSproducer));
            art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
            auto assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
                new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
            const std::vector<sim::MCShower>& inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
            const std::vector<sim::MCTrack>& inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
            auto btparts_v = common::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);

            std::vector<art::Ptr<recob::Hit>> slice_hits;
            for (const auto& pfp : slice_pfp_v) {
                auto clus_pxy_v = pfp.get<recob::Cluster>();
                for (const auto& ass_clus : clus_pxy_v) {
                    const auto& clus = clus_proxy[ass_clus.key()];
                    auto clus_hit_v = clus.get<recob::Hit>();
                    slice_hits.insert(slice_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
                }
            }
            std::sort(slice_hits.begin(), slice_hits.end());
            auto last = std::unique(slice_hits.begin(), slice_hits.end());
            slice_hits.erase(last, slice_hits.end());

            float purity = 0., completeness = 0., overlay_purity = 0.;
            int ibt = common::getAssocBtPart(slice_hits, assocMCPart, btparts_v, purity, completeness, overlay_purity);
            if (ibt >= 0) {
                slice_purity = purity;
                slice_completeness = completeness;
            }
        }
    }

    void PandoraSliceAnalysis::setBranches(TTree* _tree) {
        _tree->Branch("slice_nu_score", &slice_nu_score, "slice_nu_score/F");
        _tree->Branch("slice_purity", &slice_purity, "slice_purity/F");
        _tree->Branch("slice_completeness", &slice_completeness, "slice_completeness/F");
    }

    void PandoraSliceAnalysis::resetTTree(TTree* _tree) {
        slice_nu_score = -1.0;
        slice_purity = -1.0;
        slice_completeness = -1.0;
    }

    DEFINE_ART_CLASS_TOOL(PandoraSliceAnalysis)
}