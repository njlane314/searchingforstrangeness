#ifndef IMAGEGENERATORBASE_H
#define IMAGEGENERATORBASE_H

#include "AnalysisToolBase.h"
#include "../EventClassifier.h"
#include "../ImageProcessor.h"
#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "../CommonFunctions/Pandora.h"
#include "TVector3.h"
#include <unordered_map>
#include "../CommonFunctions/Types.h"

namespace analysis 
{
    class ImageGeneratorBase : public AnalysisToolBase {
    public:
        ImageGeneratorBase(const fhicl::ParameterSet& pset);
        virtual ~ImageGeneratorBase() = default;
        void configure(const fhicl::ParameterSet& pset) override;
        void loadBadChannelMask(const art::Event& e);
        std::pair<double, double> calculateCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits);
        std::vector<art::Ptr<recob::Hit>> collectNeutrinoHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& neutrino_slice);
        bool isNeutrinoContained(const art::Event& e, const std::vector<image::ImageProperties>& properties, const std::vector<double>& centroid_drifts);
        std::vector<image::ImageProperties> buildImageProperties(const std::vector<double>& centroid_wires, const std::vector<double>& centroid_drifts);

        virtual void analyseEvent(const art::Event& e, bool is_data) override = 0;
        virtual void analyseSlice(const art::Event& e, std::vector<common::ProxyPfpElem_t>& slicePfpVector, bool is_data, bool selected) override = 0;
        virtual void setBranches(TTree* tree) override = 0;
        virtual void resetTTree(TTree* tree) override = 0;

    protected:
        art::InputTag _MCTproducer, _WREproducer, _HITproducer, _BKTproducer, _PFPproducer, _CLSproducer,
                      _SHRproducer, _SLCproducer, _VTXproducer, _PCAproducer, _TRKproducer, _DeadChannelTag, _MCPproducer;
        std::vector<int> planes_;
        std::vector<int> width_, height_;
        float _drift_step;
        float _wire_pitch_u, _wire_pitch_v, _wire_pitch_w;
        std::vector<bool> _bad_channel_mask;
        const geo::GeometryCore* _geo;
        const detinfo::DetectorProperties* _detp;
    };

    ImageGeneratorBase::ImageGeneratorBase(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void ImageGeneratorBase::configure(const fhicl::ParameterSet& pset) {
        _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");
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
        _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
        _DeadChannelTag = pset.get<art::InputTag>("DeadChannelTag", "nfbadchannel:badchannels:OverlayDetsim");
        planes_ = pset.get<std::vector<int>>("Planes", {0, 1, 2});
        width_ = pset.get<std::vector<int>>("Width", std::vector<int>(planes_.size(), 512));
        height_ = pset.get<std::vector<int>>("Height", std::vector<int>(planes_.size(), 512));
        _geo = art::ServiceHandle<geo::Geometry>()->provider();
        _detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        auto const* det_clocks = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
        double tick_period = det_clocks->TPCClock().TickPeriod();
        double drift_velocity = _detp->DriftVelocity();
        _drift_step = tick_period * drift_velocity * 1e1;
        _wire_pitch_u = _geo->WirePitch(geo::kU);
        _wire_pitch_v = _geo->WirePitch(geo::kV);
        _wire_pitch_w = _geo->WirePitch(geo::kW);
        size_t n_channels = _geo->Nchannels();
        _bad_channel_mask.assign(n_channels + 1, false);
    }

    void ImageGeneratorBase::loadBadChannelMask(const art::Event& e) {
        art::Handle<std::vector<int>> bad_ch_h;
        if (e.getByLabel(_DeadChannelTag, bad_ch_h) && bad_ch_h.isValid()) {
            std::vector<art::Ptr<int>> bad_ch_v;
            art::fill_ptr_vector(bad_ch_v, bad_ch_h);
            for (const auto& ch : bad_ch_v) _bad_channel_mask[*ch] = true;
        }
    }

    std::pair<double, double> ImageGeneratorBase::calculateCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits) {
        double sum_charge = 0.0, sum_wire = 0.0, sum_drift = 0.0;
        for (const auto& hit : hits) {
            if (common::GetPandoraView(hit) != view) continue;
            double charge = hit->Integral();
            TVector3 hit_pos = common::GetPandoraHitPosition(e, hit, view);
            sum_charge += charge;
            sum_wire += hit_pos.Z() * charge;
            sum_drift += hit_pos.X() * charge;
        }
        return (sum_charge == 0.0) ? std::make_pair(0.0, 0.0) : std::make_pair(sum_wire / sum_charge, sum_drift / sum_charge);
    }

    std::vector<art::Ptr<recob::Hit>> ImageGeneratorBase::collectNeutrinoHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& neutrino_slice) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits;
        auto clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, _CLSproducer, proxy::withAssociated<recob::Hit>(_CLSproducer));
        for (const auto& pfp : neutrino_slice) {
            if (pfp->IsPrimary()) continue;
            for (const auto& ass_clus : pfp.get<recob::Cluster>()) {
                auto clus_hit_v = clus_proxy[ass_clus.key()].get<recob::Hit>();
                neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
            }
        }
        return neutrino_hits;
    }

    bool ImageGeneratorBase::isNeutrinoContained(const art::Event& e, const std::vector<image::ImageProperties>& properties, const std::vector<double>& centroid_drifts) {
        auto const& truth_handle = e.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
        if (truth_handle->empty() || !truth_handle->at(0).NeutrinoSet()) 
            return false;
        const auto& truth = truth_handle->at(0);
        const auto& neutrino = truth.GetNeutrino();
        const auto& neutrino_particle = neutrino.Nu();
        float corrected_vertex[3];
        common::True2RecoMappingXYZ(neutrino_particle.T(), neutrino_particle.Vx(), neutrino_particle.Vy(), neutrino_particle.Vz(), corrected_vertex);
        float vtx_u = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_U).Z();
        float vtx_v = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_V).Z();
        float vtx_w = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_W).Z();
        bool contained = true;
        for (size_t i = 0; i < planes_.size(); ++i) {
            float vtx = (planes_[i] == 0) ? vtx_u : (planes_[i] == 1) ? vtx_v : vtx_w;
            contained &= (vtx >= properties[i].origin_x() && vtx < properties[i].max_x() &&
                          centroid_drifts[i] >= properties[i].origin_y() && centroid_drifts[i] < properties[i].max_y());
        }
        return contained;
    }

    std::vector<image::ImageProperties> ImageGeneratorBase::buildImageProperties(const std::vector<double>& centroid_wires, const std::vector<double>& centroid_drifts) {
        std::vector<image::ImageProperties> properties;
        for (size_t i = 0; i < planes_.size(); ++i) {
            float pitch = (planes_[i] == 0) ? _wire_pitch_u : (planes_[i] == 1) ? _wire_pitch_v : _wire_pitch_w;
            properties.emplace_back(centroid_wires[i], centroid_drifts[i], height_[i], width_[i], pitch, _drift_step, 
                                    static_cast<geo::View_t>(planes_[i]), 0);
        }
        return properties;
    }
}

#endif