#ifndef ANALYSIS_TOPOLOGYANALYSIS_CXX
#define ANALYSIS_TOPOLOGYANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include <canvas/Persistency/Common/FindManyP.h>
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <vector>
#include <limits>
#include <map>

namespace analysis 
{
    class TopologyAnalysis : public AnalysisToolBase {
    public:
        explicit TopologyAnalysis(const fhicl::ParameterSet& pset) { configure(pset); }
        ~TopologyAnalysis() = default;

        void configure(fhicl::ParameterSet const& pset) override;
        void analyseEvent(art::Event const& e, bool fData) override;
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        art::InputTag fCLSproducer;
        art::InputTag fPFPproducer;
        art::InputTag fSLCEproducer;
        art::InputTag fHITproducer;
        int m_totaleventhits;
        double m_semanticsalience;
        int m_nhit_u, m_nhits_v, m_nhits_w;
        int m_nclusters_u, m_nclusters_v, m_nclusters_w;
        float m_charge_u, m_charge_v, m_charge_w;
        float m_wirerange_u, m_wirerange_v, m_wirerange_w;
        float m_timerange_u, m_timerange_v, m_timerange_w;

        void computeFeatures(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v);
    };

    void TopologyAnalysis::configure(fhicl::ParameterSet const& pset) {
        fCLSproducer = pset.get<art::InputTag>("CLSproducer", "pandora");
        fPFPproducer = pset.get<art::InputTag>("PFPproducer", "pandora");
        fSLCEproducer = pset.get<art::InputTag>("SLICEproducer", "pandora");
        fHITproducer = pset.get<art::InputTag>("HITproducer", "gaushit");
    }

    void TopologyAnalysis::analyseEvent(art::Event const& e, bool fData) {
        art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
        m_totaleventhits = inputHits->size();
    }

    void TopologyAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) {
        std::cout << "TopologyAnalysis::analyseSlice: " << slice_pfp_v.size() << " PFPs in slice." << std::endl;
        this->computeFeatures(e, slice_pfp_v);
    }

    void TopologyAnalysis::computeFeatures(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v) {
        if (slice_pfp_v.empty()) return;

        auto pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);

        art::FindManyP<recob::Slice> sliceAssoc(pfpHandle, e, fPFPproducer);

        size_t pfpIndex = slice_pfp_v[0].index();

        const auto& slices = sliceAssoc.at(pfpIndex);
        if (slices.empty()) return;

        const art::Ptr<recob::Slice>& slice = slices[0];

        auto sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSLCEproducer);
        art::FindManyP<recob::Hit> hitAssoc(sliceHandle, e, fSLCEproducer);
        const std::vector<art::Ptr<recob::Hit>>& sliceHits = hitAssoc.at(slice.key());

        std::map<int, std::vector<art::Ptr<recob::Hit>>> sliceHitsPerPlane;
        for (const auto& hit : sliceHits) {
            int plane = hit->WireID().planeID().Plane;
            if (plane >= 0 && plane <= 2) {
                sliceHitsPerPlane[plane].push_back(hit);
            }
        }
        int totalSliceHits = sliceHits.size();
        m_semanticsalience = (m_totaleventhits > 0) ? (static_cast<double>(totalSliceHits) * totalSliceHits) / m_totaleventhits : 0.0;
        m_nclusters_u = 0;
        m_nclusters_v = 0;
        m_nclusters_w = 0;
        for (const auto& pfp : slice_pfp_v) {
            auto clusters = pfp.get<recob::Cluster>();
            for (const auto& cluster : clusters) {
                geo::PlaneID planeID = cluster->Plane();
                if (planeID.isValid) {
                    int plane = planeID.Plane;
                    if (plane == 0) m_nclusters_u++;
                    else if (plane == 1) m_nclusters_v++;
                    else if (plane == 2) m_nclusters_w++;
                }
            }
        }
        for (int plane = 0; plane < 3; ++plane) {
            if (sliceHitsPerPlane.count(plane)) {
                const auto& hits = sliceHitsPerPlane[plane];
                int nhits = hits.size();
                float charge_sum = 0.0f;
                float min_wire = std::numeric_limits<float>::max();
                float max_wire = std::numeric_limits<float>::lowest();
                float min_time = std::numeric_limits<float>::max();
                float max_time = std::numeric_limits<float>::lowest();
                for (const auto& hit : hits) {
                    float wire = hit->WireID().Wire;
                    float time = hit->PeakTime();
                    charge_sum += hit->Integral();
                    if (wire < min_wire) min_wire = wire;
                    if (wire > max_wire) max_wire = wire;
                    if (time < min_time) min_time = time;
                    if (time > max_time) max_time = time;
                }
                float wirerange = (nhits > 0) ? max_wire - min_wire : 0.0f;
                float timerange = (nhits > 0) ? max_time - min_time : 0.0f;
                if (plane == 0) {
                    m_nhit_u = nhits;
                    m_charge_u = charge_sum;
                    m_wirerange_u = wirerange;
                    m_timerange_u = timerange;
                } else if (plane == 1) {
                    m_nhits_v = nhits;
                    m_charge_v = charge_sum;
                    m_wirerange_v = wirerange;
                    m_timerange_v = timerange;
                } else {
                    m_nhits_w = nhits;
                    m_charge_w = charge_sum;
                    m_wirerange_w = wirerange;
                    m_timerange_w = timerange;
                }
            } else {
                if (plane == 0) {
                    m_nhit_u = 0;
                    m_charge_u = 0.0f;
                    m_wirerange_u = 0.0f;
                    m_timerange_u = 0.0f;
                } else if (plane == 1) {
                    m_nhits_v = 0;
                    m_charge_v = 0.0f;
                    m_wirerange_v = 0.0f;
                    m_timerange_v = 0.0f;
                } else {
                    m_nhits_w = 0;
                    m_charge_w = 0.0f;
                    m_wirerange_w = 0.0f;
                    m_timerange_w = 0.0f;
                }
            }
        }
        std::cout << "TopologyAnalysis: U hits: " << m_nhit_u << ", V hits: " << m_nhits_v
                << ", W hits: " << m_nhits_w << std::endl;
        std::cout << "Clusters - U: " << m_nclusters_u << ", V: " << m_nclusters_v
                << ", W: " << m_nclusters_w << std::endl;
        std::cout << "Semantic Salience: " << m_semanticsalience << std::endl;
        std::cout << "Charges - U: " << m_charge_u << ", V: " << m_charge_v
                << ", W: " << m_charge_w << std::endl;
        std::cout << "Wire ranges - U: " << m_wirerange_u << ", V: " << m_wirerange_v
                << ", W: " << m_wirerange_w << std::endl;
        std::cout << "Time ranges - U: " << m_timerange_u << ", V: " << m_timerange_v
                << ", W: " << m_timerange_w << std::endl;
    }

    void TopologyAnalysis::setBranches(TTree* tree) {
        tree->Branch("nhits_u", &m_nhit_u, "nhits_u/I");
        tree->Branch("nhits_v", &m_nhits_v, "nhits_v/I");
        tree->Branch("nhits_w", &m_nhits_w, "nhits_w/I");
        tree->Branch("nclusters_u", &m_nclusters_u, "nclusters_u/I");
        tree->Branch("nclusters_v", &m_nclusters_v, "nclusters_v/I");
        tree->Branch("nclusters_w", &m_nclusters_w, "nclusters_w/I");
        tree->Branch("charge_u", &m_charge_u, "charge_u/F");
        tree->Branch("charge_v", &m_charge_v, "charge_v/F");
        tree->Branch("charge_w", &m_charge_w, "charge_w/F");
        tree->Branch("wirerange_u", &m_wirerange_u, "wirerange_u/F");
        tree->Branch("wirerange_v", &m_wirerange_v, "wirerange_v/F");
        tree->Branch("wirerange_w", &m_wirerange_w, "wirerange_w/F");
        tree->Branch("timerange_u", &m_timerange_u, "timerange_u/F");
        tree->Branch("timerange_v", &m_timerange_v, "timerange_v/F");
        tree->Branch("timerange_w", &m_timerange_w, "timerange_w/F");
        tree->Branch("semantic_salience", &m_semanticsalience, "semantic_salience/D");
        tree->Branch("totaleventhits", &m_totaleventhits, "totaleventhits/I");
    }

    void TopologyAnalysis::resetTTree(TTree* tree) {
        m_nhit_u = 0;
        m_nhits_v = 0;
        m_nhits_w = 0;
        m_nclusters_u = 0;
        m_nclusters_v = 0;
        m_nclusters_w = 0;
        m_charge_u = 0.0f;
        m_charge_v = 0.0f;
        m_charge_w = 0.0f;
        m_wirerange_u = 0.0f;
        m_wirerange_v = 0.0f;
        m_wirerange_w = 0.0f;
        m_timerange_u = 0.0f;
        m_timerange_v = 0.0f;
        m_timerange_w = 0.0f;
        m_semanticsalience = 0.0;
        m_totaleventhits = 0;
    }

    DEFINE_ART_CLASS_TOOL(TopologyAnalysis)
}

#endif