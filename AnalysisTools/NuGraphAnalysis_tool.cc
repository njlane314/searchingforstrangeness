#ifndef ANALYSIS_NUGRAPHCOUNTS_CXX
#define ANALYSIS_NUGRAPHCOUNTS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardata/Utilities/FindManyInChainP.h"

namespace analysis
{
    class NuGraphAnalysis : public AnalysisToolBase {
    public:
        explicit NuGraphAnalysis(const fhicl::ParameterSet &pset);
        virtual ~NuGraphAnalysis(){};

        void configure(fhicl::ParameterSet const &pset);

        void analyseEvent(art::Event const &e, bool fData) override;

        void analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

        void SaveTruth(art::Event const &e);

        void setBranches(TTree *_tree) override;
        
        void resetTTree(TTree *_tree) override;

    private:
        art::InputTag fCLSproducer; 
        art::InputTag fSLCproducer; 
        art::InputTag fNG2producer; 

        template <typename T, typename A>
        int arg_max(std::vector<T, A> const& vec) {
            return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
        }

        int slcng2mip;   // Number of MIP hits in the slice
        int slcng2hip;   // Number of HIP hits in the slice
        int slcng2shr;   // Number of shower hits in the slice
        int slcng2mcl;   // Number of Michel hits in the slice
        int slcng2dfs;   // Number of diffuse hits in the slice

        int clung2mip;   // Number of MIP hits in clustered hits
        int clung2hip;   // Number of HIP hits in clustered hits
        int clung2shr;   // Number of shower hits in clustered hits
        int clung2mcl;   // Number of Michel hits in clustered hits
        int clung2dfs;   // Number of diffuse hits in clustered hits

        std::vector<int> pfng2semlabel;    // Semantic label for each PFP
        std::vector<float> pfng2mipfrac;   // Fraction of MIP hits in each PFP
        std::vector<float> pfng2hipfrac;   // Fraction of HIP hits in each PFP
        std::vector<float> pfng2shrfrac;   // Fraction of shower hits in each PFP
        std::vector<float> pfng2mclfrac;   // Fraction of Michel hits in each PFP
        std::vector<float> pfng2dfsfrac;   // Fraction of diffuse hits in each PFP
        std::vector<float> pfng2mipavrg;   // Average MIP score for each PFP 
        std::vector<float> pfng2hipavrg;   // Average HIP score for each PFP
        std::vector<float> pfng2shravrg;   // Average shower score for each PFP
        std::vector<float> pfng2mclavrg;   // Average Michel score for each PFP
        std::vector<float> pfng2dfsavrg;   // Average diffuse score for each PFP
    };

    NuGraphAnalysis::NuGraphAnalysis(const fhicl::ParameterSet &p) {
        this->configure(p);
    }

    void NuGraphAnalysis::configure(fhicl::ParameterSet const &p) {
        fCLSproducer = p.get<art::InputTag>("CLSproducer");
        fSLCproducer = p.get<art::InputTag>("SLCproducer");
        fNG2producer = p.get<art::InputTag>("NG2producer");
    }

    void NuGraphAnalysis::analyseEvent(art::Event const &e, bool fData) {}

    void NuGraphAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) {
        common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
        art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
        auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));

        // auto GNNDescription = e.getHandle<anab::MVADescription<5>>(art::InputTag("NuGraph", "semantic"));
        auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(e, art::InputTag("gaushit"), 
                                    //proxy::withParallelData<anab::FeatureVector<1>>(art::InputTag("NuGraph", "filter")),
                                    proxy::withParallelData<anab::FeatureVector<5>>(fNG2producer));

        std::vector<int> ng2semclucounts(5,0);
        std::vector<int> ng2semslccounts(5,0);
        for (auto& h : hitsWithScores) {
            auto scores = h.get<anab::FeatureVector<5>>();
            std::vector<float> ng2semscores;
            for (size_t i=0;i<scores.size();i++) ng2semscores.push_back(scores[i]);
            unsigned int sem_label = arg_max(ng2semscores);
            ng2semslccounts[sem_label]++;
        }

        for (auto pfp : slice_pfp_v) {
            if (pfp->PdgCode()!=11 && pfp->PdgCode()!=13) continue;

            std::vector<art::Ptr<recob::Hit>> hit_v;
            auto clus_pxy_v = pfp.get<recob::Cluster>();
            if (clus_pxy_v.size() != 0) {
                for (auto ass_clus : clus_pxy_v) {
                    const auto &clus = clus_proxy[ass_clus.key()];
                    auto clus_hit_v = clus.get<recob::Hit>();
                    for (const auto &hit : clus_hit_v)
                    hit_v.push_back(hit);
                } 
            }

            if (hit_v.size()>0) {
                std::vector<int> ng2sempfpcounts(5,0);
                std::vector<float> ng2sempfptotscr(5,0);
                for (auto& hit : hit_v) {
                    auto scores = hitsWithScores[hit.key()].get<anab::FeatureVector<5>>();
                    std::vector<float> ng2semscores;
                    for (size_t i=0;i<scores.size();i++) {
                        ng2semscores.push_back(scores[i]);
                        ng2sempfptotscr[i] += scores[i];
                    }
                    unsigned int sem_label = arg_max(ng2semscores);
                    ng2sempfpcounts[sem_label]++;
                    ng2semclucounts[sem_label]++;
                }
                pfng2semlabel.push_back(arg_max(ng2sempfpcounts));
                pfng2mipfrac.push_back(float(ng2sempfpcounts[0])/hit_v.size());
                pfng2hipfrac.push_back(float(ng2sempfpcounts[1])/hit_v.size());
                pfng2shrfrac.push_back(float(ng2sempfpcounts[2])/hit_v.size());
                pfng2mclfrac.push_back(float(ng2sempfpcounts[3])/hit_v.size());
                pfng2dfsfrac.push_back(float(ng2sempfpcounts[4])/hit_v.size());
                pfng2mipavrg.push_back(float(ng2sempfptotscr[0])/hit_v.size());
                pfng2hipavrg.push_back(float(ng2sempfptotscr[1])/hit_v.size());
                pfng2shravrg.push_back(float(ng2sempfptotscr[2])/hit_v.size());
                pfng2mclavrg.push_back(float(ng2sempfptotscr[3])/hit_v.size());
                pfng2dfsavrg.push_back(float(ng2sempfptotscr[4])/hit_v.size());
            } else {
                pfng2semlabel.push_back(-1);
                pfng2mipfrac.push_back(-1);
                pfng2hipfrac.push_back(-1);
                pfng2shrfrac.push_back(-1);
                pfng2mclfrac.push_back(-1);
                pfng2dfsfrac.push_back(-1);
                pfng2mipavrg.push_back(-1);
                pfng2hipavrg.push_back(-1);
                pfng2shravrg.push_back(-1);
                pfng2mclavrg.push_back(-1);
                pfng2dfsavrg.push_back(-1);
            }
        }   
        
        slcng2mip = ng2semslccounts[0];
        slcng2hip = ng2semslccounts[1];
        slcng2shr = ng2semslccounts[2];
        slcng2mcl = ng2semslccounts[3];
        slcng2dfs = ng2semslccounts[4];
        
        clung2mip = ng2semclucounts[0];
        clung2hip = ng2semclucounts[1];
        clung2shr = ng2semclucounts[2];
        clung2mcl = ng2semclucounts[3];
        clung2dfs = ng2semclucounts[4];

        return;
    }

    void NuGraphAnalysis::setBranches(TTree *_tree) {
        _tree->Branch("slcng2mip", &slcng2mip, "slcng2mip/I");
        _tree->Branch("slcng2hip", &slcng2hip, "slcng2hip/I");
        _tree->Branch("slcng2shr", &slcng2shr, "slcng2shr/I");
        _tree->Branch("slcng2mcl", &slcng2mcl, "slcng2mcl/I");
        _tree->Branch("slcng2dfs", &slcng2dfs, "slcng2dfs/I");
        
        _tree->Branch("clung2mip", &clung2mip, "clung2mip/I");
        _tree->Branch("clung2hip", &clung2hip, "clung2hip/I");
        _tree->Branch("clung2shr", &clung2shr, "clung2shr/I");
        _tree->Branch("clung2mcl", &clung2mcl, "clung2mcl/I");
        _tree->Branch("clung2dfs", &clung2dfs, "clung2dfs/I");
        
        _tree->Branch("pfng2semlabel", &pfng2semlabel);
        _tree->Branch("pfng2mipfrac", &pfng2mipfrac);
        _tree->Branch("pfng2hipfrac", &pfng2hipfrac);
        _tree->Branch("pfng2shrfrac", &pfng2shrfrac);
        _tree->Branch("pfng2mclfrac", &pfng2mclfrac);
        _tree->Branch("pfng2dfsfrac", &pfng2dfsfrac);
        _tree->Branch("pfng2mipavrg", &pfng2mipavrg);
        _tree->Branch("pfng2hipavrg", &pfng2hipavrg);
        _tree->Branch("pfng2shravrg", &pfng2shravrg);
        _tree->Branch("pfng2mclavrg", &pfng2mclavrg);
        _tree->Branch("pfng2dfsavrg", &pfng2dfsavrg);
    }

    void NuGraphAnalysis::resetTTree(TTree *_tree) {
        slcng2mip = std::numeric_limits<int>::min();
        slcng2hip = std::numeric_limits<int>::min();
        slcng2shr = std::numeric_limits<int>::min();
        slcng2mcl = std::numeric_limits<int>::min();
        slcng2dfs = std::numeric_limits<int>::min();
      
        clung2mip = std::numeric_limits<int>::min();
        clung2hip = std::numeric_limits<int>::min();
        clung2shr = std::numeric_limits<int>::min();
        clung2mcl = std::numeric_limits<int>::min();
        clung2dfs = std::numeric_limits<int>::min();
       
        pfng2semlabel.clear();
        pfng2mipfrac.clear();
        pfng2hipfrac.clear();
        pfng2shrfrac.clear();
        pfng2mclfrac.clear();
        pfng2dfsfrac.clear();
        pfng2mipavrg.clear();
        pfng2hipavrg.clear();
        pfng2shravrg.clear();
        pfng2mclavrg.clear();
        pfng2dfsavrg.clear();
    }

    DEFINE_ART_CLASS_TOOL(NuGraphAnalysis)
} 

#endif