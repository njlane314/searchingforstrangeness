#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "../CommonDefs/Image.h"
#include "../CommonDefs/Pandora.h"
#include "../CommonDefs/Types.h"
#include "../CommonDefs/ReconstructionLabels.h"
#include "../CommonDefs/PrimaryLabels.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>
#include "canvas/Persistency/Common/FindManyP.h"
#include <vector>
#include <string>
#include <map>
#include "TTree.h"
#include <cmath>

namespace analysis 
{
    class ImageAnalysis : public AnalysisToolBase {
    public:
        explicit ImageAnalysis(fhicl::ParameterSet const& p);
        virtual ~ImageAnalysis() = default;

        void configure(const fhicl::ParameterSet& p) override;

        void analyseEvent(art::Event const& e, bool _is_data) override {}

        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) override;

        void setBranches(TTree* _tree) override {}

        void resetTTree(TTree* _tree) override {}

    private:
        art::InputTag fPFPproducer;
        art::InputTag fCLSproducer;
        art::InputTag fHITproducer;
        art::InputTag fWIREproducer;
        art::InputTag fMCPproducer;
        std::string fBackTrackerLabel;
        double fGammaThreshold;
        double fHadronThreshold;

        int _image_width;
        int _image_height;
        float _adc_image_threshold;

        const geo::GeometryCore* _geo;
        const detinfo::DetectorProperties* _detp;

        float _drift_step;
        float _wire_pitch_u;
        float _wire_pitch_v;
        float _wire_pitch_w;

        art::TFileService* _tfs;
        TTree* _image_tree;

        int nrows;
        int ncols;

        std::vector<std::vector<float>> raw_images;
        std::vector<std::vector<int>> nugraph_images;
        std::vector<std::vector<int>> truth_images;
        bool is_vertex_image; 

        std::vector<art::Ptr<recob::Hit>> collectNeutrinoHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);
        std::pair<double, double> calculateCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits);
        void constructImages(const art::Event& e, const std::vector<image::ImageProperties>& properties, 
                             std::vector<image::Image<float>>& raw_images, 
                             std::vector<image::Image<int>>& nugraph_images, 
                             std::vector<image::Image<int>>& truth_images);
    };

    ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void ImageAnalysis::configure(const fhicl::ParameterSet& p) {
        fPFPproducer = p.get<art::InputTag>("PFPproducer");
        fCLSproducer = p.get<art::InputTag>("CLSproducer");
        fHITproducer = p.get<art::InputTag>("HITproducer");
        fWIREproducer = p.get<art::InputTag>("WIREproducer");
        fMCPproducer = p.get<art::InputTag>("MCPproducer");
        fBackTrackerLabel = p.get<std::string>("BackTrackerLabel", "gaushit");
        fGammaThreshold = p.get<double>("GammaThreshold", 0.1);
        fHadronThreshold = p.get<double>("HadronThreshold", 0.1);
        _image_width = p.get<int>("ImageWidth", 512);
        _image_height = p.get<int>("ImageHeight", 512);
        _adc_image_threshold = p.get<float>("ADCthreshold", 4.0);

        _geo = art::ServiceHandle<geo::Geometry>()->provider();
        _detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();

        double tick_period = clock->TPCClock().TickPeriod();
        double drift_velocity = _detp->DriftVelocity();
        _drift_step = tick_period * drift_velocity * 1e1;
        _wire_pitch_u = _geo->WirePitch(geo::kU);
        _wire_pitch_v = _geo->WirePitch(geo::kV);
        _wire_pitch_w = _geo->WirePitch(geo::kW);
        _tfs = &(*art::ServiceHandle<art::TFileService>());

        _image_tree = _tfs->make<TTree>("ImageTree", "Tree containing image data");

        _image_tree->Branch("nrows", &nrows);
        _image_tree->Branch("ncols", &ncols);
        _image_tree->Branch("raw_images", &raw_images);
        _image_tree->Branch("nugraph_images", &nugraph_images);
        _image_tree->Branch("truth_images", &truth_images);
        _image_tree->Branch("is_vertex_image", &is_vertex_image, "is_vertex_image/O");

        raw_images.resize(3);
        nugraph_images.resize(3);
        truth_images.resize(3);
    }

    void ImageAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits = this->collectNeutrinoHits(e, slice_pfp_v);
        auto [centroid_wire_u, centroid_drift_u] = this->calculateCentroid(e, common::TPC_VIEW_U, neutrino_hits);
        auto [centroid_wire_v, centroid_drift_v] = this->calculateCentroid(e, common::TPC_VIEW_V, neutrino_hits);
        auto [centroid_wire_w, centroid_drift_w] = this->calculateCentroid(e, common::TPC_VIEW_W, neutrino_hits);
        std::vector<image::ImageProperties> properties;
        properties.emplace_back(centroid_wire_u, centroid_drift_u, _image_height, _image_width, _wire_pitch_u, _drift_step, geo::kU);
        properties.emplace_back(centroid_wire_v, centroid_drift_v, _image_height, _image_width, _wire_pitch_v, _drift_step, geo::kV);
        properties.emplace_back(centroid_wire_w, centroid_drift_w, _image_height, _image_width, _wire_pitch_w, _drift_step, geo::kW);

        double nu_vtx[3] = {};
        bool has_neutrino = false;
        for (const auto& pfp : slice_pfp_v) {
            if (pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 12 || std::abs(pfp->PdgCode()) == 14)) {
                auto vtx_v = pfp.get<recob::Vertex>();
                if (vtx_v.size() == 1) {
                    vtx_v[0]->XYZ(nu_vtx);
                    has_neutrino = true;
                    break;
                }
            }
        }
        if (!has_neutrino) {
            is_vertex_image = false;
        } else {
            is_vertex_image = true;
            for (size_t view_id = 0; view_id < 3; ++view_id) {
                geo::View_t view = properties[view_id].view();
                double wire_coord;
                if (view == geo::kW) {
                    wire_coord = nu_vtx[2]; // z
                } else if (view == geo::kU) {
                    wire_coord = nu_vtx[2] * std::cos(1.04719758034) - nu_vtx[1] * std::sin(1.04719758034);
                } else if (view == geo::kV) {
                    wire_coord = nu_vtx[2] * std::cos(-1.04719758034) - nu_vtx[1] * std::sin(-1.04719758034);
                } else {
                    continue;
                }
                double x = nu_vtx[0];
                size_t row = properties[view_id].row(x);
                size_t col = properties[view_id].col(wire_coord);
                if (row >= properties[view_id].height() || col >= properties[view_id].width()) {
                    is_vertex_image = false;
                    break;
                }
            }
        }

        std::vector<image::Image<float>> raw_img_vec;
        std::vector<image::Image<int>> nugraph_img_vec;
        std::vector<image::Image<int>> truth_img_vec;
        this->constructImages(e, properties, raw_img_vec, nugraph_img_vec, truth_img_vec);

        nrows = properties[0].height();
        ncols = properties[0].width();

        raw_images.clear();
        nugraph_images.clear();
        truth_images.clear();

        for (size_t view_id = 0; view_id < 3; ++view_id) {
            raw_images.push_back(raw_img_vec[view_id].data());
            nugraph_images.push_back(nugraph_img_vec[view_id].data());
            truth_images.push_back(truth_img_vec[view_id].data());
        }

        _image_tree->Fill();
    }

    std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectNeutrinoHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits;
        auto clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
        for (const auto& pfp : pfp_pxy_v) {
            if (pfp->IsPrimary()) continue;
            for (auto ass_clus : pfp.get<recob::Cluster>()) {
                auto clus_hit_v = clus_proxy[ass_clus.key()].get<recob::Hit>();
                neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
            }
        }
        return neutrino_hits;
    }

    std::pair<double, double> ImageAnalysis::calculateCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits) {
        double sum_charge = 0.0;
        double sum_wire = 0.0;
        double sum_drift = 0.0;
        for (const auto& hit : hits) {
            if (common::GetPandoraView(hit) != view) continue;
            double charge = hit->Integral();
            TVector3 hit_pos = common::GetPandoraHitPosition(e, hit, view);
            sum_charge += charge;
            sum_wire += hit_pos.Z() * charge;
            sum_drift += hit_pos.X() * charge;
        }
        if (sum_charge == 0.0) return {0.0, 0.0};
        return {sum_wire / sum_charge, sum_drift / sum_charge};
    }

    void ImageAnalysis::constructImages(const art::Event& e,
                                            const std::vector<image::ImageProperties>& properties,
                                            std::vector<image::Image<float>>& raw_images,
                                            std::vector<image::Image<int>>& nugraph_images,
                                            std::vector<image::Image<int>>& truth_images) {
        raw_images.clear();
        nugraph_images.clear();
        truth_images.clear();

        for (const auto& prop : properties) {
            image::Image<float> raw_image(prop);
            raw_image.clear(0.0);
            raw_images.push_back(std::move(raw_image));

            image::Image<int> nugraph_image(prop);
            nugraph_image.clear(static_cast<int>(reco_labels::ReconstructionLabel::empty));
            nugraph_images.push_back(std::move(nugraph_image));

            image::Image<int> truth_image(prop);
            truth_image.clear(static_cast<int>(truth_labels::PrimaryLabel::empty));
            truth_images.push_back(std::move(truth_image));
        }

        auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(fWIREproducer);
        auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);

        art::Handle<std::vector<simb::MCParticle>> mcpHandle;
        bool hasMCParticles = e.getByLabel(fMCPproducer, mcpHandle);

        art::FindManyP<recob::Hit> wire_hit_assoc(wireHandle, e, fHITproducer);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hitHandle, e, fBackTrackerLabel);

        std::vector<truth_labels::PrimaryLabel> primary_labels;
        std::vector<reco_labels::ReconstructionLabel> reco_labels;

        if (hasMCParticles && mcpHandle.isValid()) {
            primary_labels = truth_labels::classifyParticles(e, fMCPproducer, fGammaThreshold, fHadronThreshold);
            reco_labels = reco_labels::classifyParticles(e, fMCPproducer, fGammaThreshold, fHadronThreshold);
        } else {
            primary_labels = {truth_labels::PrimaryLabel::cosmic};
            reco_labels = {reco_labels::ReconstructionLabel::cosmic};
        }

        std::map<int, size_t> trackid_to_index;
        if (hasMCParticles && mcpHandle.isValid()) {
            for (size_t i = 0; i < mcpHandle->size(); ++i) {
                trackid_to_index[mcpHandle->at(i).TrackId()] = i;
            }
        }

        for (size_t wire_idx = 0; wire_idx < wireHandle->size(); ++wire_idx) {
            const auto& wire = wireHandle->at(wire_idx);
            auto ch_id = wire.Channel();
            std::vector<geo::WireID> wire_ids = _geo->ChannelToWire(ch_id);
            if (wire_ids.empty()) continue;
            geo::View_t view = _geo->View(wire_ids.front().planeID());
            size_t view_idx = view - geo::kU;

            const geo::WireGeo* wire_geo = _geo->WirePtr(wire_ids.front());
            TVector3 center = wire_geo->GetCenter();
            TVector3 wire_center(center.X(), center.Y(), center.Z());
            double wire_coord = (view == geo::kW) ? wire_center.Z() :
                                (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                                    (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));

            auto hits = wire_hit_assoc.at(wire_idx);

            std::sort(hits.begin(), hits.end(), [](const art::Ptr<recob::Hit>& a, const art::Ptr<recob::Hit>& b) {
                return a->StartTick() < b->StartTick();
            });

            size_t hit_idx = 0;

            for (const auto& range : wire.SignalROI().get_ranges()) {
                const auto& adcs = range.data();
                int start_tick = range.begin_index();
                for (size_t idx = 0; idx < adcs.size(); ++idx) {
                    int tick = start_tick + idx;
                    double x = _detp->ConvertTicksToX(tick, wire_ids.front().planeID());
                    size_t row = properties[view_idx].row(x);
                    size_t col = properties[view_idx].col(wire_coord);
                    if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                    reco_labels::ReconstructionLabel reco_label;
                    truth_labels::PrimaryLabel primary_label;

                    if (!hasMCParticles) {
                        reco_label = reco_labels::ReconstructionLabel::cosmic;
                        primary_label = truth_labels::PrimaryLabel::cosmic;
                    } else {
                        reco_label = reco_labels::ReconstructionLabel::invisible;
                        primary_label = truth_labels::PrimaryLabel::other;

                        while (hit_idx < hits.size() && hits[hit_idx]->EndTick() <= tick) {
                            ++hit_idx;
                        }
                        if (hit_idx < hits.size() && hits[hit_idx]->StartTick() <= tick && tick < hits[hit_idx]->EndTick()) {
                            const auto& hit = hits[hit_idx];
                            auto bkth_data = mcp_bkth_assoc.data(hit.key());
                            if (!bkth_data.empty()) {
                                for (size_t i = 0; i < bkth_data.size(); ++i) {
                                    if (bkth_data[i]->isMaxIDE == 1) {
                                        int track_id = mcp_bkth_assoc.at(hit.key())[i]->TrackId();
                                        auto it = trackid_to_index.find(track_id);
                                        if (it != trackid_to_index.end()) {
                                            size_t particle_idx = it->second;
                                            if (particle_idx < reco_labels.size()) {
                                                reco_label = reco_labels[particle_idx];
                                            }
                                            if (particle_idx < primary_labels.size()) {
                                                primary_label = primary_labels[particle_idx];
                                            }
                                        }
                                        break;
                                    }
                                }
                            } else {
                                reco_label = reco_labels::ReconstructionLabel::cosmic;
                                primary_label = truth_labels::PrimaryLabel::cosmic;
                            }
                        }
                    }

                    if (adcs[idx] > _adc_image_threshold) {
                        raw_images[view_idx].set(row, col, adcs[idx]);
                        nugraph_images[view_idx].set(row, col, static_cast<int>(reco_label), false);
                        truth_images[view_idx].set(row, col, static_cast<int>(primary_label), false);
                    }
                }
            }
        }
    }

    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}