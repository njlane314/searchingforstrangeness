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
#include <fstream>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TVector3.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/System/TriggerNamesService.h"

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
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/GeometryUtilities.h"  
#include "lardata/Utilities/FindManyInChainP.h"

#include "AnalysisToolBase.h"
#include "../CommonDefs/Types.h"
#include "../CommonDefs/Pandora.h"
#include "../ImageAlgorithm.h"

#ifdef ClassDef
#undef ClassDef
#endif
#include "torch/torch.h"
#include "torch/script.h"

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
        art::InputTag fWIREproducer;
        art::InputTag fHITproducer;
        art::InputTag fBKTproducer;
        art::InputTag fMCPproducer;
        art::InputTag fCLSproducer; 

        std::unique_ptr<ImageAlgorithm> image_algo_;
        std::vector<Image<float>> detector_image_tensor_;
        std::vector<Image<int>> semantic_image_tensor_;

        void initialiseImageProperties(const art::Event& event, std::vector<ImageProperties>& image_props, const std::vector<art::Ptr<recob::Hit>>& slice_hits);

        std::vector<float> _detector_plane_image_u;
        std::vector<float> _detector_plane_image_v;
        std::vector<float> _detector_plane_image_w;
        std::vector<int> _semantic_plane_image_u;
        std::vector<int> _semantic_plane_image_v;
        std::vector<int> _semantic_plane_image_w;
    };

    ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void ImageAnalysis::configure(const fhicl::ParameterSet& pset) {
        image_algo_ = std::make_unique<ImageAlgorithm>(pset);

        fWIREproducer = pset.get<art::InputTag>("WireProducer", "butcher");
        fHITproducer = pset.get<art::InputTag>("HitProducer", "gaushit");
        fBKTproducer = pset.get<art::InputTag>("BackTrackerProducer", "gaushitTruthMatch");
        fMCPproducer = pset.get<art::InputTag>("MCParticleProducer", "largeant");
        fCLSproducer = pset.get<art::InputTag>("ClusterProducer", "pandora");  
    }

    void ImageAnalysis::setBranches(TTree* _tree) {
        _tree->Branch("detector_image_u", &_detector_plane_image_u);
        _tree->Branch("detector_image_v", &_detector_plane_image_v);
        _tree->Branch("detector_image_w", &_detector_plane_image_w);
        _tree->Branch("semantic_image_u", &_semantic_plane_image_u);
        _tree->Branch("semantic_image_v", &_semantic_plane_image_v);
        _tree->Branch("semantic_image_w", &_semantic_plane_image_w);
    }

    void ImageAnalysis::resetTTree(TTree* /*tree*/) {
        _detector_plane_image_u.clear();
        _detector_plane_image_v.clear();
        _detector_plane_image_w.clear();
        _semantic_plane_image_u.clear();
        _semantic_plane_image_v.clear();
        _semantic_plane_image_w.clear();
    }

    void ImageAnalysis::initialiseImageProperties(
        const art::Event& event,
        std::vector<ImageProperties>& image_props,
        const std::vector<art::Ptr<recob::Hit>>& slice_hits) {
        auto geo = art::ServiceHandle<geo::Geometry>()->provider();
        size_t n_planes = geo->Nplanes();

        image_props.clear();
        image_props.reserve(n_planes);

        auto calculateCentroid = [&](geo::View_t view, const std::vector<art::Ptr<recob::Hit>>& hits) {
            double sum_charge = 0.0;
            double weighted_sum_drift = 0.0;
            double weighted_sum_wire = 0.0;

            for (const auto& hit : hits) {
                if (hit->View() != view) continue;
                const double charge = hit->Integral();
                auto pos = common::GetPandoraHitPosition(event, hit, common::GetPandoraView(hit));
                double drift = pos.X();
                double wire = pos.Z();
                weighted_sum_drift += drift * charge;
                weighted_sum_wire += wire * charge;
                sum_charge += charge;
            }

            if (sum_charge > 0.0) {
                double centroid_drift = weighted_sum_drift / sum_charge;
                double centroid_wire = weighted_sum_wire / sum_charge;
                std::cout << "Centroid for view " << view << ": Drift = " << centroid_drift
                          << ", Wire = " << centroid_wire << std::endl;
                return std::make_pair(centroid_drift, centroid_wire);
            } else {
                return std::make_pair(0.0, 0.0);
            }
        };

        for (size_t i = 0; i < n_planes; ++i) {
            geo::View_t view = static_cast<geo::View_t>(i);
            auto [centroid_drift, centroid_wire] = calculateCentroid(view, slice_hits);
            if (centroid_drift == 0.0 && centroid_wire == 0.0) continue;

            size_t width = image_algo_->getImageWidth();
            size_t height = image_algo_->getImageHeight();
            double pixel_w = image_algo_->getDriftStep();
            double pixel_h = image_algo_->getWirePitch(view);

            double center_x = centroid_drift;
            double center_y = centroid_wire;

            image_props.emplace_back(center_x, center_y, width, height, pixel_h, pixel_w, view);
        }
    }

    void ImageAnalysis::analyseSlice(
        art::Event const& event,
        std::vector<common::ProxyPfpElem_t>& slice_pfp_v,
        bool is_data,
        bool selected) {
        std::vector<art::Ptr<recob::Wire>> wires;
        if (auto wireHandle = event.getValidHandle<std::vector<recob::Wire>>(fWIREproducer))
            art::fill_ptr_vector(wires, wireHandle);

        std::vector<art::Ptr<recob::Hit>> hits;
        if (auto hitHandle = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer))
            art::fill_ptr_vector(hits, hitHandle);

        std::vector<art::Ptr<simb::MCParticle>> mcps;
        if (auto mcpHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer))
            art::fill_ptr_vector(mcps, mcpHandle);

        std::vector<art::Ptr<recob::Hit>> slice_hits;
        auto clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(event, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
        for (const auto& pfp : slice_pfp_v) {
            if (pfp->IsPrimary()) continue;
            for (auto ass_clus : pfp.get<recob::Cluster>()) {
                auto clus_hit_v = clus_proxy[ass_clus.key()].get<recob::Hit>();
                slice_hits.insert(slice_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
            }
        }

        std::vector<ImageProperties> image_props;
        detector_image_tensor_.clear();
        semantic_image_tensor_.clear();
        this->initialiseImageProperties(event, image_props, slice_hits);
        if (image_algo_ && !image_props.empty()) {
            image_algo_->generateViewImages(event, image_props, wires, mcps);
            detector_image_tensor_ = image_algo_->getDetectorImages();
            semantic_image_tensor_ = image_algo_->getSemanticImages();
        }

        for (size_t i = 0; i < detector_image_tensor_.size() && i < image_props.size(); ++i) {
            const auto& image = detector_image_tensor_[i];
            const auto& semantic_image = semantic_image_tensor_[i];
            const auto& prop = image_props[i];
            if (prop.view() == geo::kU) {
                _detector_plane_image_u = image.data();
                _semantic_plane_image_u = semantic_image.data();
            } else if (prop.view() == geo::kV) {
                _detector_plane_image_v = image.data();
                _semantic_plane_image_v = semantic_image.data();
            } else if (prop.view() == geo::kW || prop.view() == geo::kY) {
                _detector_plane_image_w = image.data();
                _semantic_plane_image_w = semantic_image.data();
            }
        }
    }

    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif