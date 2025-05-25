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

#include "../Image.h"
#include "../ImageClassifiers.h"

namespace analysis 
{
    class ImageAnalysis : public AnalysisToolBase {
    public:
        explicit ImageAnalysis(fhicl::ParameterSet const& p);
        virtual ~ImageAnalysis() = default;

        void configure(const fhicl::ParameterSet& p) override;
        void analyseEvent(art::Event const& e, bool _is_data) override {}
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) override;
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;

    private:
        art::InputTag fPFPproducer;
        art::InputTag fCLSproducer;
        art::InputTag fHITproducer;
        art::InputTag fWIREproducer;
        art::InputTag fMCPproducer;
        std::string fBackTrackerLabel;

        bool fProcessMC;
        bool fUseCheatRecoLabels;

        double fGammaThreshold;
        double fHadronThreshold;
        double fLeptonThreshold;

        int _image_width;
        int _image_height;
        float _adc_image_threshold;

        const geo::GeometryCore* _geo;
        const detinfo::DetectorProperties* _detp;

        float _drift_step;
        float _wire_pitch_u;
        float _wire_pitch_v;
        float _wire_pitch_w;

        TTree* _tree = nullptr;

        std::vector<float> _raw_image_u;
        std::vector<float> _raw_image_v;
        std::vector<float> _raw_image_w;
        std::vector<int> _reco_image_u;
        std::vector<int> _reco_image_v;
        std::vector<int> _reco_image_w;
        std::vector<int> _true_image_u;
        std::vector<int> _true_image_v;
        std::vector<int> _true_image_w;

        std::unique_ptr<RecoLabelClassifier> fRecoClassifier;
        std::unique_ptr<TruthLabelClassifier> fTruthClassifier;

        std::vector<art::Ptr<recob::Hit>> collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);
        std::pair<double, double> calculateChargeCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits);
        void constructPixelImages(const art::Event& e, const std::vector<ImageProperties>& properties,
                                std::vector<Image<float>>& raw_images,
                                std::vector<Image<int>>& reco_images,
                                std::vector<Image<int>>& true_images);
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
        fLeptonThreshold = p.get<double>("LeptonThreshold", 0.1);

        fProcessMC = p.get<bool>("ProcessMC", true);
        fUseCheatRecoLabels = p.get<bool>("UseCheatRecoLabels", true);

        _image_width = p.get<int>("ImageWidth", 512);
        _image_height = p.get<int>("ImageHeight", 512);
        _adc_image_threshold = p.get<float>("ADCthreshold", 1.0);

        _geo = art::ServiceHandle<geo::Geometry>()->provider();
        _detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();

        double tick_period = clock->TPCClock().TickPeriod();
        double drift_velocity = _detp->DriftVelocity();
        _drift_step = tick_period * drift_velocity * 1e1;
        _wire_pitch_u = _geo->WirePitch(geo::kU);
        _wire_pitch_v = _geo->WirePitch(geo::kV);
        _wire_pitch_w = _geo->WirePitch(geo::kW);

        if (fProcessMC) {
            fTruthClassifier = std::make_unique<TruthLabelClassifier>(fMCPproducer);
            fRecoClassifier = std::make_unique<RecoLabelClassifier>(fMCPproducer, fGammaThreshold, fHadronThreshold, fUseCheatRecoLabels);
        }
    }

    void ImageAnalysis::setBranches(TTree* _tree_ptr) {
        _tree = _tree_ptr;
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

    void ImageAnalysis::resetTTree(TTree* _tree_ptr) {
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
        std::cout << "image slice" << std::endl;
        std::vector<art::Ptr<recob::Hit>> neutrino_hits = this->collectSliceHits(e, slice_pfp_v);
        auto [centroid_wire_u, centroid_drift_u] = this->calculateChargeCentroid(e, common::TPC_VIEW_U, neutrino_hits);
        auto [centroid_wire_v, centroid_drift_v] = this->calculateChargeCentroid(e, common::TPC_VIEW_V, neutrino_hits);
        auto [centroid_wire_w, centroid_drift_w] = this->calculateChargeCentroid(e, common::TPC_VIEW_W, neutrino_hits);
        std::vector<ImageProperties> properties;
        properties.emplace_back(centroid_wire_u, centroid_drift_u, _image_height, _image_width, _drift_step, _wire_pitch_u, geo::kU);
        properties.emplace_back(centroid_wire_v, centroid_drift_v, _image_height, _image_width, _drift_step, _wire_pitch_v, geo::kV);
        properties.emplace_back(centroid_wire_w, centroid_drift_w, _image_height, _image_width, _drift_step, _wire_pitch_w, geo::kW);

        std::vector<Image<float>> out_raw_images;
        std::vector<Image<int>> out_reco_images;
        std::vector<Image<int>> out_true_images;
        this->constructPixelImages(e, properties, out_raw_images, out_reco_images, out_true_images);

        _raw_image_u.clear();
        _raw_image_v.clear();
        _raw_image_w.clear();
        _reco_image_u.clear();
        _reco_image_v.clear();
        _reco_image_w.clear();
        _true_image_u.clear();
        _true_image_v.clear();
        _true_image_w.clear();

        _raw_image_u = out_raw_images[0].data();  // U view (geo::kU)
        _raw_image_v = out_raw_images[1].data();  // V view (geo::kV)
        _raw_image_w = out_raw_images[2].data();  // W view (geo::kW)
        _reco_image_u = out_reco_images[0].data();
        _reco_image_v = out_reco_images[1].data();
        _reco_image_w = out_reco_images[2].data();
        _true_image_u = out_true_images[0].data();
        _true_image_v = out_true_images[1].data();
        _true_image_w = out_true_images[2].data();

        if (_tree) _tree->Fill();

        std::cout << "finished image slice" << std::endl;
    }

    std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
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

    std::pair<double, double> ImageAnalysis::calculateChargeCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits) {
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

    void ImageAnalysis::constructPixelImages(const art::Event& e,
                                            const std::vector<ImageProperties>& properties,
                                            std::vector<Image<float>>& out_raw_images,
                                            std::vector<Image<int>>& out_reco_images,
                                            std::vector<Image<int>>& out_true_images) {
        out_raw_images.clear();
        out_reco_images.clear();
        out_true_images.clear();

        for (const auto& prop : properties) {
            Image<float> raw_image(prop);
            raw_image.clear(0.0);
            out_raw_images.push_back(std::move(raw_image));

            Image<int> reco_image(prop);
            reco_image.clear(static_cast<int>(RecoLabelClassifier::ReconstructionLabel::Empty));
            out_reco_images.push_back(std::move(reco_image));

            Image<int> true_image(prop);
            true_image.clear(static_cast<int>(TruthLabelClassifier::TruthPrimaryLabel::Empty));
            out_true_images.push_back(std::move(true_image));
        }

        auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(fWIREproducer);
        auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);

        art::Handle<std::vector<simb::MCParticle>> mcpHandle;
        bool hasMCParticles = e.getByLabel(fMCPproducer, mcpHandle);

        art::FindManyP<recob::Hit> wire_hit_assoc(wireHandle, e, fHITproducer);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hitHandle, e, fBackTrackerLabel);

        std::vector<TruthLabelClassifier::TruthPrimaryLabel> classified_true_labels;
        std::vector<RecoLabelClassifier::ReconstructionLabel> classified_reco_labels;

        if (fProcessMC && hasMCParticles && mcpHandle.isValid() && fTruthClassifier && fRecoClassifier) {
            classified_true_labels = fTruthClassifier->classifyParticles(e);
            classified_reco_labels = fRecoClassifier->classifyParticles(e);
        }

        std::map<int, size_t> trackid_to_index;
        if (fProcessMC && hasMCParticles && mcpHandle.isValid()) {
            for (size_t i = 0; i < mcpHandle->size(); ++i) {
                trackid_to_index[mcpHandle->at(i).TrackId()] = i;
            }
        }

        // Calibration factors for U, V, W planes
        std::vector<double> calibration_factors = {232.0, 249.0, 243.7};

        for (size_t wire_idx = 0; wire_idx < wireHandle->size(); ++wire_idx) {
            const auto& wire = wireHandle->at(wire_idx);
            auto ch_id = wire.Channel();
            std::vector<geo::WireID> wire_ids = _geo->ChannelToWire(ch_id);
            if (wire_ids.empty()) continue;
            geo::View_t view = _geo->View(wire_ids.front().planeID());
            size_t view_idx = static_cast<size_t>(view); 

            const geo::WireGeo* wire_geo = _geo->WirePtr(wire_ids.front());
            TVector3 center = wire_geo->GetCenter();
            TVector3 wire_center(center.X(), center.Y(), center.Z());
            double wire_coord = (view == geo::kW) ? wire_center.Z() :
                                (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));

            auto hits_for_wire = wire_hit_assoc.at(wire_idx);
            std::vector<art::Ptr<recob::Hit>> sorted_hits = hits_for_wire;

            std::sort(sorted_hits.begin(), sorted_hits.end(), [](const art::Ptr<recob::Hit>& a, const art::Ptr<recob::Hit>& b) {
                return a->StartTick() < b->StartTick();
            });

            size_t current_hit_idx_in_sorted_list = 0;

            for (const auto& range : wire.SignalROI().get_ranges()) {
                const auto& adcs = range.data();
                int start_tick = range.begin_index();
                for (size_t idx_in_adc_range = 0; idx_in_adc_range < adcs.size(); ++idx_in_adc_range) {
                    int tick = start_tick + idx_in_adc_range;
                    double x_drift_pos = _detp->ConvertTicksToX(static_cast<double>(tick), wire_ids.front().planeID());
                    size_t row = properties[view_idx].row(x_drift_pos);
                    size_t col = properties[view_idx].col(wire_coord);
                    if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                    RecoLabelClassifier::ReconstructionLabel current_reco_label_for_pixel;
                    TruthLabelClassifier::TruthPrimaryLabel current_true_label_for_pixel;

                    if (fProcessMC && hasMCParticles && mcpHandle.isValid() && fTruthClassifier && fRecoClassifier) {
                        current_reco_label_for_pixel = RecoLabelClassifier::ReconstructionLabel::Invisible;
                        current_true_label_for_pixel = TruthLabelClassifier::TruthPrimaryLabel::Other;

                        while (current_hit_idx_in_sorted_list < sorted_hits.size() && sorted_hits[current_hit_idx_in_sorted_list]->EndTick() <= tick) {
                            ++current_hit_idx_in_sorted_list;
                        }

                        if (current_hit_idx_in_sorted_list < sorted_hits.size() &&
                            sorted_hits[current_hit_idx_in_sorted_list]->StartTick() <= tick &&
                            tick < sorted_hits[current_hit_idx_in_sorted_list]->EndTick()) {
                            const art::Ptr<recob::Hit>& matched_hit = sorted_hits[current_hit_idx_in_sorted_list];
                            
                            std::vector<art::Ptr<simb::MCParticle>> mcp_particles_ass_to_hit;
                            std::vector<anab::BackTrackerHitMatchingData const*> bkth_data_ass_to_hit;
                            mcp_bkth_assoc.get(matched_hit.key(), mcp_particles_ass_to_hit, bkth_data_ass_to_hit);

                            if (!bkth_data_ass_to_hit.empty()) {
                                for (size_t i_bkth = 0; i_bkth < bkth_data_ass_to_hit.size(); ++i_bkth) {
                                    if (bkth_data_ass_to_hit[i_bkth] && bkth_data_ass_to_hit[i_bkth]->isMaxIDE) {
                                        int track_id = mcp_particles_ass_to_hit[i_bkth]->TrackId();
                                        auto it_trackid = trackid_to_index.find(track_id);
                                        if (it_trackid != trackid_to_index.end()) {
                                            size_t particle_mcp_idx = it_trackid->second;
                                            if (particle_mcp_idx < classified_reco_labels.size()) {
                                                current_reco_label_for_pixel = classified_reco_labels[particle_mcp_idx];
                                            }
                                            if (particle_mcp_idx < classified_true_labels.size()) {
                                                current_true_label_for_pixel = classified_true_labels[particle_mcp_idx];
                                            }
                                        }
                                        break; 
                                    }
                                }
                            } else {
                                current_reco_label_for_pixel = RecoLabelClassifier::ReconstructionLabel::Cosmic;
                                current_true_label_for_pixel = TruthLabelClassifier::TruthPrimaryLabel::Cosmic;
                            }
                        }
                    } else {
                        current_reco_label_for_pixel = RecoLabelClassifier::ReconstructionLabel::Cosmic;
                        current_true_label_for_pixel = TruthLabelClassifier::TruthPrimaryLabel::Cosmic;
                    }

                    if (adcs[idx_in_adc_range] > _adc_image_threshold) {
                        float adc = adcs[idx_in_adc_range];
                        double energy = static_cast<double>(adc) * calibration_factors[view_idx];
                        double log_energy = std::log10(energy);
                        out_raw_images[view_idx].set(row, col, static_cast<float>(log_energy));
                        out_reco_images[view_idx].set(row, col, static_cast<int>(current_reco_label_for_pixel), false);
                        out_true_images[view_idx].set(row, col, static_cast<int>(current_true_label_for_pixel), false);
                    }
                }
            }
        }
    }

    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif