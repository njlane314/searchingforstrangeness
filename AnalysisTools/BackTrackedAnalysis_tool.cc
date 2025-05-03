#ifndef ANALYSIS_BACKTRACKEDANALYSIS_CXX
#define ANALYSIS_BACKTRACKEDANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/ProximityClustering.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/CalibrationFuncs.h"

namespace analysis {
    class BackTrackedAnalysis : public AnalysisToolBase {

    public:
        BackTrackedAnalysis(const fhicl::ParameterSet &pset);
        ~BackTrackedAnalysis() {};
        void configure(fhicl::ParameterSet const &pset);
        void analyseEvent(art::Event const &e, bool fData) override;
        void analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;
        void setBranches(TTree *_tree) override;
        void resetTTree(TTree *_tree) override;

    private:
        art::InputTag fCLSproducer;
        art::InputTag fMCRproducer;
        art::InputTag fBacktrackTag;
        art::InputTag fHproducer;

        std::vector<int> bt_pfp_pdg;
        std::vector<int> bt_pfp_nhits;
        std::vector<int> bt_pfp_nplanehits_U;
        std::vector<int> bt_pfp_nplanehits_V;
        std::vector<int> bt_pfp_nplanehits_Y;
        std::vector<int> bt_pfp_planesubclusters_U;
        std::vector<int> bt_pfp_planesubclusters_V;
        std::vector<int> bt_pfp_planesubclusters_Y;
        std::vector<float> bt_pfp_planesubhitfracmax_U;
        std::vector<float> bt_pfp_planesubhitfracmax_V;
        std::vector<float> bt_pfp_planesubhitfracmax_Y;
        std::vector<float> bt_pfp_vtx_x;
        std::vector<float> bt_pfp_vtx_y;
        std::vector<float> bt_pfp_vtx_z;

        std::vector<int> bt_backtracked_pdg;
        std::vector<float> bt_backtracked_e;
        std::vector<int> bt_backtracked_tid;
        std::vector<float> bt_backtracked_purity;
        std::vector<float> bt_backtracked_completeness;
        std::vector<float> bt_backtracked_overlay_purity;
        std::vector<float> bt_backtracked_px;
        std::vector<float> bt_backtracked_py;
        std::vector<float> bt_backtracked_pz;
        std::vector<float> bt_backtracked_start_x;
        std::vector<float> bt_backtracked_start_y;
        std::vector<float> bt_backtracked_start_z;
        std::vector<float> bt_backtracked_start_t;
        std::vector<float> bt_backtracked_start_U;
        std::vector<float> bt_backtracked_start_V;
        std::vector<float> bt_backtracked_start_Y;
        std::vector<float> bt_backtracked_sce_start_x;
        std::vector<float> bt_backtracked_sce_start_y;
        std::vector<float> bt_backtracked_sce_start_z;
        std::vector<float> bt_backtracked_sce_start_U;
        std::vector<float> bt_backtracked_sce_start_V;
        std::vector<float> bt_backtracked_sce_start_Y;
    };

    BackTrackedAnalysis::BackTrackedAnalysis(const fhicl::ParameterSet &p) {
        fCLSproducer = p.get<art::InputTag>("CLSproducer");
        fMCRproducer = p.get<art::InputTag>("MCRproducer");
        fBacktrackTag = p.get<art::InputTag>("BacktrackTag");
        fHproducer = p.get<art::InputTag>("Hproducer");
    }

    void BackTrackedAnalysis::configure(fhicl::ParameterSet const &p) {
    }

    void BackTrackedAnalysis::analyseEvent(art::Event const &e, bool fData) {
    }

    void BackTrackedAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) {
        common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                              proxy::withAssociated<recob::Hit>(fCLSproducer));

        std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
        std::vector<common::BtPart> btparts_v;
        if (!fData) {
            const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
            const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
            art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
            assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
                new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
            btparts_v = common::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);
        }

        for (auto pfp : slice_pfp_v) {
            if (pfp->IsPrimary()) continue;

            bt_pfp_pdg.push_back(pfp->PdgCode());

            auto vertices = pfp.get<recob::Vertex>();
            if (vertices.size() == 1) {
                double xyz[3];
                vertices.at(0)->XYZ(xyz);
                bt_pfp_vtx_x.push_back(xyz[0]);
                bt_pfp_vtx_y.push_back(xyz[1]);
                bt_pfp_vtx_z.push_back(xyz[2]);
            } else {
                bt_pfp_vtx_x.push_back(std::numeric_limits<float>::lowest());
                bt_pfp_vtx_y.push_back(std::numeric_limits<float>::lowest());
                bt_pfp_vtx_z.push_back(std::numeric_limits<float>::lowest());
            }

            std::vector<art::Ptr<recob::Hit>> hit_v;
            auto clus_pxy_v = pfp.get<recob::Cluster>();

            int nplanehits_U = 0, nplanehits_V = 0, nplanehits_Y = 0;
            int planesubclusters_U = 0, planesubclusters_V = 0, planesubclusters_Y = 0;
            float planesubhitfracmax_U = 0, planesubhitfracmax_V = 0, planesubhitfracmax_Y = 0;

            for (auto ass_clus : clus_pxy_v) {
                const auto &clus = clus_proxy[ass_clus.key()];
                auto clus_hit_v = clus.get<recob::Hit>();
                auto nhits = clus_hit_v.size();
                std::vector<art::Ptr<recob::Hit>> cluster_hits_v;
                for (size_t h = 0; h < clus_hit_v.size(); h++) {
                    cluster_hits_v.push_back(clus_hit_v[h]);
                }
                int nclus = 0;
                float hitfracmax = 0.;
                std::vector<std::vector<unsigned int>> out_cluster_v;
                if (nhits) {
                    common::cluster(cluster_hits_v, out_cluster_v, 2.0, 1.0);
                    for (size_t nc = 0; nc < out_cluster_v.size(); nc++) {
                        auto clus_hit_idx_v = out_cluster_v.at(nc);
                        int nhitclus = clus_hit_idx_v.size();
                        if (nhitclus > 3.)
                            nclus += 1;
                        float hitfrac = nhitclus / nhits;
                        if (hitfrac > hitfracmax)
                            hitfracmax = hitfrac;
                    }
                }
                if (clus->Plane().Plane == 0) {
                    nplanehits_U += nhits;
                    planesubclusters_U += nclus;
                    planesubhitfracmax_U = hitfracmax;
                } else if (clus->Plane().Plane == 1) {
                    nplanehits_V += nhits;
                    planesubclusters_V += nclus;
                    planesubhitfracmax_V = hitfracmax;
                } else if (clus->Plane().Plane == 2) {
                    nplanehits_Y += nhits;
                    planesubclusters_Y += nclus;
                    planesubhitfracmax_Y = hitfracmax;
                }
                for (const auto &hit : clus_hit_v) {
                    hit_v.push_back(hit);
                }
            }
            bt_pfp_nhits.push_back(hit_v.size());
            bt_pfp_nplanehits_U.push_back(nplanehits_U);
            bt_pfp_nplanehits_V.push_back(nplanehits_V);
            bt_pfp_nplanehits_Y.push_back(nplanehits_Y);
            bt_pfp_planesubclusters_U.push_back(planesubclusters_U);
            bt_pfp_planesubclusters_V.push_back(planesubclusters_V);
            bt_pfp_planesubclusters_Y.push_back(planesubclusters_Y);
            bt_pfp_planesubhitfracmax_U.push_back(planesubhitfracmax_U);
            bt_pfp_planesubhitfracmax_V.push_back(planesubhitfracmax_V);
            bt_pfp_planesubhitfracmax_Y.push_back(planesubhitfracmax_Y);

            if (!fData) {
                float purity = 0., completeness = 0., overlay_purity = 0.;
                int ibt = common::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness, overlay_purity);
                if (ibt >= 0) {
                    auto &mcp = btparts_v[ibt];
                    auto PDG = mcp.pdg;
                    bt_backtracked_e.push_back(mcp.e);
                    bt_backtracked_tid.push_back(mcp.tids.at(0));
                    bt_backtracked_pdg.push_back(PDG);
                    bt_backtracked_purity.push_back(purity);
                    bt_backtracked_completeness.push_back(completeness);
                    bt_backtracked_overlay_purity.push_back(overlay_purity);
                    bt_backtracked_px.push_back(mcp.px);
                    bt_backtracked_py.push_back(mcp.py);
                    bt_backtracked_pz.push_back(mcp.pz);
                    bt_backtracked_start_x.push_back(mcp.start_x);
                    bt_backtracked_start_y.push_back(mcp.start_y);
                    bt_backtracked_start_z.push_back(mcp.start_z);
                    bt_backtracked_start_t.push_back(mcp.start_t);
                    bt_backtracked_start_U.push_back(common::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 0));
                    bt_backtracked_start_V.push_back(common::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 1));
                    bt_backtracked_start_Y.push_back(common::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 2));
                    float reco_st[3] = {mcp.start_x, mcp.start_y, mcp.start_z};
                    if (PDG == 11 || PDG == 22) {
                        reco_st[0] += common::x_offset(mcp.start_t);
                    } else {
                        common::True2RecoMappingXYZ(mcp.start_t, mcp.start_x, mcp.start_y, mcp.start_z, reco_st);
                    }
                    bt_backtracked_sce_start_x.push_back(reco_st[0]);
                    bt_backtracked_sce_start_y.push_back(reco_st[1]);
                    bt_backtracked_sce_start_z.push_back(reco_st[2]);
                    bt_backtracked_sce_start_U.push_back(common::YZtoPlanecoordinate(reco_st[1], reco_st[2], 0));
                    bt_backtracked_sce_start_V.push_back(common::YZtoPlanecoordinate(reco_st[1], reco_st[2], 1));
                    bt_backtracked_sce_start_Y.push_back(common::YZtoPlanecoordinate(reco_st[1], reco_st[2], 2));
                } else {
                    bt_backtracked_e.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_tid.push_back(std::numeric_limits<int>::lowest());
                    bt_backtracked_pdg.push_back(std::numeric_limits<int>::lowest());
                    bt_backtracked_purity.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_completeness.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_overlay_purity.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_px.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_py.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_pz.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_start_x.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_start_y.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_start_z.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_start_t.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_start_U.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_start_V.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_start_Y.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_sce_start_x.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_sce_start_y.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_sce_start_z.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_sce_start_U.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_sce_start_V.push_back(std::numeric_limits<float>::lowest());
                    bt_backtracked_sce_start_Y.push_back(std::numeric_limits<float>::lowest());
                }
            } else {
                bt_backtracked_e.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_tid.push_back(std::numeric_limits<int>::lowest());
                bt_backtracked_pdg.push_back(std::numeric_limits<int>::lowest());
                bt_backtracked_purity.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_completeness.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_overlay_purity.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_px.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_py.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_pz.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_start_x.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_start_y.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_start_z.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_start_t.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_start_U.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_start_V.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_start_Y.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_sce_start_x.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_sce_start_y.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_sce_start_z.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_sce_start_U.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_sce_start_V.push_back(std::numeric_limits<float>::lowest());
                bt_backtracked_sce_start_Y.push_back(std::numeric_limits<float>::lowest());
            }
        }
    }

    void BackTrackedAnalysis::setBranches(TTree *_tree) {
        _tree->Branch("bt_pfp_pdg", "std::vector<int>", &bt_pfp_pdg);
        _tree->Branch("bt_pfp_nhits", "std::vector<int>", &bt_pfp_nhits);
        _tree->Branch("bt_pfp_nplanehits_U", "std::vector<int>", &bt_pfp_nplanehits_U);
        _tree->Branch("bt_pfp_nplanehits_V", "std::vector<int>", &bt_pfp_nplanehits_V);
        _tree->Branch("bt_pfp_nplanehits_Y", "std::vector<int>", &bt_pfp_nplanehits_Y);
        _tree->Branch("bt_pfp_planesubclusters_U", "std::vector<int>", &bt_pfp_planesubclusters_U);
        _tree->Branch("bt_pfp_planesubclusters_V", "std::vector<int>", &bt_pfp_planesubclusters_V);
        _tree->Branch("bt_pfp_planesubclusters_Y", "std::vector<int>", &bt_pfp_planesubclusters_Y);
        _tree->Branch("bt_pfp_planesubhitfracmax_U", "std::vector<float>", &bt_pfp_planesubhitfracmax_U);
        _tree->Branch("bt_pfp_planesubhitfracmax_V", "std::vector<float>", &bt_pfp_planesubhitfracmax_V);
        _tree->Branch("bt_pfp_planesubhitfracmax_Y", "std::vector<float>", &bt_pfp_planesubhitfracmax_Y);
        _tree->Branch("bt_pfp_vtx_x", "std::vector<float>", &bt_pfp_vtx_x);
        _tree->Branch("bt_pfp_vtx_y", "std::vector<float>", &bt_pfp_vtx_y);
        _tree->Branch("bt_pfp_vtx_z", "std::vector<float>", &bt_pfp_vtx_z);
        _tree->Branch("bt_backtracked_pdg", "std::vector<int>", &bt_backtracked_pdg);
        _tree->Branch("bt_backtracked_e", "std::vector<float>", &bt_backtracked_e);
        _tree->Branch("bt_backtracked_tid", "std::vector<int>", &bt_backtracked_tid);
        _tree->Branch("bt_backtracked_purity", "std::vector<float>", &bt_backtracked_purity);
        _tree->Branch("bt_backtracked_completeness", "std::vector<float>", &bt_backtracked_completeness);
        _tree->Branch("bt_backtracked_overlay_purity", "std::vector<float>", &bt_backtracked_overlay_purity);
        _tree->Branch("bt_backtracked_px", "std::vector<float>", &bt_backtracked_px);
        _tree->Branch("bt_backtracked_py", "std::vector<float>", &bt_backtracked_py);
        _tree->Branch("bt_backtracked_pz", "std::vector<float>", &bt_backtracked_pz);
        _tree->Branch("bt_backtracked_start_x", "std::vector<float>", &bt_backtracked_start_x);
        _tree->Branch("bt_backtracked_start_y", "std::vector<float>", &bt_backtracked_start_y);
        _tree->Branch("bt_backtracked_start_z", "std::vector<float>", &bt_backtracked_start_z);
        _tree->Branch("bt_backtracked_start_t", "std::vector<float>", &bt_backtracked_start_t);
        _tree->Branch("bt_backtracked_start_U", "std::vector<float>", &bt_backtracked_start_U);
        _tree->Branch("bt_backtracked_start_V", "std::vector<float>", &bt_backtracked_start_V);
        _tree->Branch("bt_backtracked_start_Y", "std::vector<float>", &bt_backtracked_start_Y);
        _tree->Branch("bt_backtracked_sce_start_x", "std::vector<float>", &bt_backtracked_sce_start_x);
        _tree->Branch("bt_backtracked_sce_start_y", "std::vector<float>", &bt_backtracked_sce_start_y);
        _tree->Branch("bt_backtracked_sce_start_z", "std::vector<float>", &bt_backtracked_sce_start_z);
        _tree->Branch("bt_backtracked_sce_start_U", "std::vector<float>", &bt_backtracked_sce_start_U);
        _tree->Branch("bt_backtracked_sce_start_V", "std::vector<float>", &bt_backtracked_sce_start_V);
        _tree->Branch("bt_backtracked_sce_start_Y", "std::vector<float>", &bt_backtracked_sce_start_Y);
    }

    void BackTrackedAnalysis::resetTTree(TTree *_tree) {
        bt_pfp_pdg.clear();
        bt_pfp_nhits.clear();
        bt_pfp_nplanehits_U.clear();
        bt_pfp_nplanehits_V.clear();
        bt_pfp_nplanehits_Y.clear();
        bt_pfp_planesubclusters_U.clear();
        bt_pfp_planesubclusters_V.clear();
        bt_pfp_planesubclusters_Y.clear();
        bt_pfp_planesubhitfracmax_U.clear();
        bt_pfp_planesubhitfracmax_V.clear();
        bt_pfp_planesubhitfracmax_Y.clear();
        bt_pfp_vtx_x.clear();
        bt_pfp_vtx_y.clear();
        bt_pfp_vtx_z.clear();
        bt_backtracked_e.clear();
        bt_backtracked_tid.clear();
        bt_backtracked_pdg.clear();
        bt_backtracked_purity.clear();
        bt_backtracked_completeness.clear();
        bt_backtracked_overlay_purity.clear();
        bt_backtracked_px.clear();
        bt_backtracked_py.clear();
        bt_backtracked_pz.clear();
        bt_backtracked_start_x.clear();
        bt_backtracked_start_y.clear();
        bt_backtracked_start_z.clear();
        bt_backtracked_start_t.clear();
        bt_backtracked_start_U.clear();
        bt_backtracked_start_V.clear();
        bt_backtracked_start_Y.clear();
        bt_backtracked_sce_start_x.clear();
        bt_backtracked_sce_start_y.clear();
        bt_backtracked_sce_start_z.clear();
        bt_backtracked_sce_start_U.clear();
        bt_backtracked_sce_start_V.clear();
        bt_backtracked_sce_start_Y.clear();
    }

    DEFINE_ART_CLASS_TOOL(BackTrackedAnalysis)
} // namespace analysis

#endif