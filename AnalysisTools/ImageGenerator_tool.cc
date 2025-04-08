#ifndef ANALYSIS_IMAGEGENERATOR_CXX
#define ANALYSIS_IMAGEGENERATOR_CXX

#include <torch/torch.h>
#include <torch/script.h>
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Region.h"
#include "CommonFunctions/Types.h"
#include "CommonFunctions/Geometry.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "SignatureTools/SignatureToolBase.h"
#include "ImageProcessor.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include "ImageProcessor.h"
#include "AnalysisToolBase.h"

namespace analysis
{
    class ImageGenerator : public AnalysisToolBase {
    public:
        explicit ImageGenerator(const fhicl::ParameterSet& pset);
        ~ImageGenerator() = default;

        ImageGenerator(const ImageGenerator&) = delete;
        ImageGenerator(ImageGenerator&&) = delete;
        ImageGenerator& operator=(const ImageGenerator&) = delete;
        ImageGenerator& operator=(ImageGenerator&&) = delete;

        void configure(const fhicl::ParameterSet& pset);
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;
        void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slicePfpVector, bool is_data, bool selected) override;
        void analyseEvent(const art::Event& event, bool is_data) override;

    private:
        art::InputTag _MCTproducer;
        art::InputTag _MCPproducer;

        std::vector<int> planes_;
        std::vector<int> width_; 
        std::vector<int> height_;
        std::vector<std::vector<float>> input_;
        std::vector<float> energy_fraction_in_image_;
        bool is_neutrino_in_image_;

        int _image_width, _image_height;
        float _drift_step;
        float _wire_pitch_u, _wire_pitch_v, _wire_pitch_w;
        float vtx_u_, vtx_v_, vtx_w_;

        std::vector<bool> _bad_channel_mask;
        const geo::GeometryCore* _geo;
        const detinfo::DetectorProperties* _detp;
        art::InputTag _WREproducer, _HITproducer, _BKTproducer, _PFPproducer, _CLSproducer, _SHRproducer, _SLCproducer, _VTXproducer, _PCAproducer, _TRKproducer, _DeadChannelTag;

        void buildPFPMap(const common::ProxyPfpColl_t& pfp_pxy_col);
        void addDaughters(const common::ProxyPfpElem_t& pfp_pxy, const common::ProxyPfpColl_t& pfp_pxy_col, std::vector<common::ProxyPfpElem_t>& slice_v);
        std::vector<art::Ptr<recob::Hit>> collectNeutrinoHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& neutrino_slice);
        std::pair<double, double> calculateCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits);
        bool isVertexInImage(const std::vector<image::ImageProperties>& properties, float vtx_u, float vtx_v, float vtx_w);
    };

    ImageGenerator::ImageGenerator(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void ImageGenerator::configure(const fhicl::ParameterSet& pset) {
        _image_width = pset.get<int>("ImageWidth", 512);
        _image_height = pset.get<int>("ImageHeight", 512);
        _DeadChannelTag = pset.get<art::InputTag>("DeadChannelTag", "nfbadchannel:badchannels:OverlayDetsim");
        _WREproducer = pset.get<art::InputTag>("WREproducer", "butcher");
        _HITproducer = pset.get<art::InputTag>("HITproducer", "gaushit");
        _BKTproducer = pset.get<art::InputTag>("BKTproducer", "gaushitTruthMatch");
        _PFPproducer = pset.get<art::InputTag>("PFPproducer", "pandora");
        _CLSproducer = pset.get<art::InputTag>("CLSproducer", "pandora");
        _SHRproducer = pset.get<art::InputTag>("SHRproducer", "pandora");
        _SLCproducer = pset.get<art::InputTag>("SLCproducer", "pandora");
        _VTXproducer = pset.get<art::InputTag>("VTXproducer", "pandora");
        _PCAproducer = pset.get<art::InputTag>("PCAproducer", "pandora");
        _TRKproducer = pset.get<art::InputTag>("TRKproducer", "pandora");

        _geo = art::ServiceHandle<geo::Geometry>()->provider();
        auto const* det_clocks = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
        _detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        double tick_period = det_clocks->TPCClock().TickPeriod();
        double drift_velocity = _detp->DriftVelocity();
        _drift_step = tick_period * drift_velocity * 1e1; 
        _wire_pitch_u = _geo->WirePitch(geo::kU);
        _wire_pitch_v = _geo->WirePitch(geo::kV);
        _wire_pitch_w = _geo->WirePitch(geo::kW);
    }

    void ImageGenerator::analyseEvent(art::Event const& e, bool is_data) {
        if (is_data)
            return;

        auto const& truth_handle = e.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);   
        if (truth_handle->empty() || !truth_handle->at(0).NeutrinoSet()) 
            return; 
        
        const auto& truth = truth_handle->at(0);
        const auto& neutrino = truth.GetNeutrino();
        const auto& neutrino_particle = neutrino.Nu();
        float corrected_vertex[3];
        common::True2RecoMappingXYZ(neutrino_particle.T(), neutrino_particle.Vx(), neutrino_particle.Vy(), neutrino_particle.Vz(), corrected_vertex);
        vtx_u_ = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_U).Z();
        vtx_v_ = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_V).Z();
        vtx_w_ = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_W).Z();
    }

    void ImageGenerator::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) {
        if (is_data) 
            return; 

        size_t n_channels = _geo->Nchannels();
        _bad_channel_mask.assign(n_channels + 1, false);
        art::Handle<std::vector<int>> bad_ch_h;
        if (e.getByLabel(_DeadChannelTag, bad_ch_h) && bad_ch_h.isValid()) {
            std::vector<art::Ptr<int>> bad_ch_v;
            art::fill_ptr_vector(bad_ch_v, bad_ch_h);
            for (const auto& ch : bad_ch_v) {
                _bad_channel_mask[*ch] = true;
            }
        }

        std::vector<art::Ptr<recob::Wire>> wire_vec;
        if (auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(_WREproducer)) 
            art::fill_ptr_vector(wire_vec, wireHandle);

        auto neutrino_hits = this->collectNeutrinoHits(e, slice_pfp_v);
        auto [centroid_wire_u, centroid_drift_u] = this->calculateCentroid(e, common::TPC_VIEW_U, neutrino_hits);
        auto [centroid_wire_v, centroid_drift_v] = this->calculateCentroid(e, common::TPC_VIEW_V, neutrino_hits);
        auto [centroid_wire_w, centroid_drift_w] = this->calculateCentroid(e, common::TPC_VIEW_W, neutrino_hits);

        std::vector<image::ImageProperties> properties;
        properties.emplace_back(centroid_wire_u, centroid_drift_u, _image_height, _image_width, _wire_pitch_u, _drift_step, geo::kU, 0);
        properties.emplace_back(centroid_wire_v, centroid_drift_v, _image_height, _image_width, _wire_pitch_v, _drift_step, geo::kV, 0);
        properties.emplace_back(centroid_wire_w, centroid_drift_w, _image_height, _image_width, _wire_pitch_w, _drift_step, geo::kW, 0);

        input_ = image::extractImages(image::constructInputImages(properties, wire_vec, *_geo, _bad_channel_mask));

        std::vector<double> total_energy_per_view(3, 0.0);
        std::vector<double> energy_in_image_per_view(3, 0.0); 
        auto const* det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        for (const auto& hit : neutrino_hits) {
            geo::View_t view = hit->View(); 
            size_t view_idx = (view == geo::kU) ? 0 : (view == geo::kV) ? 1 : 2;
            double energy = hit->Integral();
            total_energy_per_view[view_idx] += energy;
            double x = det_props->ConvertTicksToX(hit->PeakTime(), hit->WireID()); 

            TVector3 wire_center = _geo->Cryostat(hit->WireID().Cryostat)
                                    .TPC(hit->WireID().TPC)
                                    .Plane(hit->WireID().Plane)
                                    .Wire(hit->WireID().Wire)
                                    .GetCenter();
            float coord = (view == geo::kW) ? wire_center.Z() :
                      (view == geo::kU) ? (wire_center.Z() * cos(1.04719758034) - wire_center.Y() * sin(1.04719758034)) :
                                          (wire_center.Z() * cos(-1.04719758034) - wire_center.Y() * sin(-1.04719758034));

            const auto& prop = properties[view_idx];
            size_t row = prop.row(x);
            size_t col = prop.col(coord);
            if (row != static_cast<size_t>(-1) && col != static_cast<size_t>(-1)) {
                energy_in_image_per_view[view_idx] += energy;
            }
        }
        for (size_t i = 0; i < 3; ++i) {
            energy_fraction_in_image_[i] = (total_energy_per_view[i] > 0.0) ?
                static_cast<float>(energy_in_image_per_view[i] / total_energy_per_view[i]) : 0.0f;
        }

        bool in_image_u = (vtx_u_ >= properties[0].origin_x() && vtx_u_ < properties[0].max_x() &&
                           centroid_drift_u >= properties[0].origin_y() && centroid_drift_u < properties[0].max_y());
        bool in_image_v = (vtx_v_ >= properties[1].origin_x() && vtx_v_ < properties[1].max_x() &&
                           centroid_drift_v >= properties[1].origin_y() && centroid_drift_v < properties[1].max_y());
        bool in_image_w = (vtx_w_ >= properties[2].origin_x() && vtx_w_ < properties[2].max_x() &&
                           centroid_drift_w >= properties[2].origin_y() && centroid_drift_w < properties[2].max_y());
        is_neutrino_in_image_ = in_image_u && in_image_v && in_image_w;
    }

    std::vector<art::Ptr<recob::Hit>> ImageGenerator::collectNeutrinoHits(const art::Event& e,
                                                                        const std::vector<common::ProxyPfpElem_t>& neutrino_slice) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits;
        auto clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(
            e, _CLSproducer, proxy::withAssociated<recob::Hit>(_CLSproducer));
        for (const auto& pfp : neutrino_slice) {
            if (pfp->IsPrimary()) continue; 
            for (const auto& ass_clus : pfp.get<recob::Cluster>()) {
                auto clus_hit_v = clus_proxy[ass_clus.key()].get<recob::Hit>();
                neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
            }
        }
        return neutrino_hits;
    }

    std::pair<double, double> ImageGenerator::calculateCentroid(const art::Event& e, common::PandoraView view,
                                                            const std::vector<art::Ptr<recob::Hit>>& hits) {
        double sum_charge = 0.0, sum_wire = 0.0, sum_drift = 0.0;
        for (const auto& hit : hits) {
            if (common::GetPandoraView(hit) != view) continue;
            double charge = hit->Integral();
            TVector3 hit_pos = common::GetPandoraHitPosition(e, hit, view);
            sum_charge += charge;
            sum_wire += hit_pos.Z() * charge;  
            sum_drift += hit_pos.X() * charge; 
        }
        return (sum_charge == 0.0) ? std::make_pair(0.0, 0.0) : 
            std::make_pair(sum_wire / sum_charge, sum_drift / sum_charge);
    }

    void ImageGenerator::setBranches(TTree* tree) {
        tree->Branch("input", &input_);
        tree->Branch("input_contained_frac", &energy_fraction_in_image_);
        tree->Branch("is_neutrino_in_image_", &is_neutrino_in_image_);
    }

    void ImageGenerator::resetTTree(TTree* /*tree*/) {
        input_.clear();
        energy_fraction_in_image_.assign(3, 0.0f);
        is_neutrino_in_image_ = false;
        vtx_u_ = std::numeric_limits<float>::lowest();
        vtx_v_ = std::numeric_limits<float>::lowest();
        vtx_w_ = std::numeric_limits<float>::lowest();
    }

    DEFINE_ART_CLASS_TOOL(ImageGenerator)
}

#endif