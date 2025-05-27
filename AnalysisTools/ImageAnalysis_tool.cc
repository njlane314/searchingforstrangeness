#ifndef ANALYSIS_IMAGE_CXX
#define ANALYSIS_IMAGE_CXX

#include <vector>
#include <string>
#include <map>
#include <array>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <limits>
#include <optional>

#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TVector3.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>
#include "lardataobj/AnalysisBase/MVAOutput.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "AnalysisToolBase.h"
#include "../CommonDefs/Pandora.h"
#include "../CommonDefs/Types.h"

#include "../ImageAlgorithm.h"

namespace analysis 
{
    class ImageAnalysis : public AnalysisToolBase {
    public:
        explicit ImageAnalysis(fhicl::ParameterSet const& p);
        virtual ~ImageAnalysis() = default;

        void configure(const fhicl::ParameterSet& p) override;
        void analyseEvent(art::Event const& e, bool _is_data) override {}
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        art::InputTag fCLSproducer;

        std::vector<float> _raw_image_u;
        std::vector<float> _raw_image_v;
        std::vector<float> _raw_image_w;
        std::vector<int> _reco_image_u;
        std::vector<int> _reco_image_v;
        std::vector<int> _reco_image_w;
        std::vector<int> _true_image_u;
        std::vector<int> _true_image_v;
        std::vector<int> _true_image_w;

        std::unique_ptr<ImageAlgorithm> image_algo_;
    };

    ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void ImageAnalysis::configure(const fhicl::ParameterSet& p) {
        fCLSproducer = p.get<art::InputTag>("CLSproducer");

        image_algo_ = std::make_unique<ImageAlgorithm>(p);
    }

    void ImageAnalysis::setBranches(TTree* _tree) {
        _tree->Branch("raw_image_u", &_raw_image_u);
        _tree->Branch("raw_image_v", &_raw_image_v);
        _tree->Branch("raw_image_w", &_raw_image_w);
        _tree->Branch("reco_image_u", &_reco_image_u);
        _tree->Branch("reco_image_v", &_reco_image_v);
        _tree->Branch("reco_image_w", &_reco_image_w);
        _tree->Branch("true_image_u", &_true_image_u);
        _tree->Branch("true_image_v", &_true_image_v);
        _tree->Branch("true_image_w", &_true_image_w);
    }

    void ImageAnalysis::resetTTree(TTree* _tree) {
        _raw_image_u.clear();
        _raw_image_v.clear();
        _raw_image_w.clear();
        _reco_image_u.clear();
        _reco_image_v.clear();
        _reco_image_w.clear();
        _true_image_u.clear();
        _true_image_v.clear();
        _true_image_w.clear();
    }

    void ImageAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits;
        auto clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
        for (const auto& pfp : slice_pfp_v) {
            if (pfp->IsPrimary()) continue;
            for (auto ass_clus : pfp.get<recob::Cluster>()) {
                auto clus_hit_v = clus_proxy[ass_clus.key()].get<recob::Hit>();
                neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
            }
        }

        std::vector<Image<float>> out_raw_images;
        std::vector<Image<int>> out_reco_images;
        std::vector<Image<int>> out_true_images;
        image_algo_->generateSliceImages(e, neutrino_hits, out_raw_images, out_reco_images, out_true_images);

        _raw_image_u.clear();
        _raw_image_v.clear();
        _raw_image_w.clear();
        _reco_image_u.clear();
        _reco_image_v.clear();
        _reco_image_w.clear();
        _true_image_u.clear();
        _true_image_v.clear();
        _true_image_w.clear();

        _raw_image_u = out_raw_images[0].data();
        _raw_image_v = out_raw_images[1].data();
        _raw_image_w = out_raw_images[2].data();
        _reco_image_u = out_reco_images[0].data();
        _reco_image_v = out_reco_images[1].data();
        _reco_image_w = out_reco_images[2].data();
        _true_image_u = out_true_images[0].data();
        _true_image_v = out_true_images[1].data();
        _true_image_w = out_true_images[2].data();
    }

    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif