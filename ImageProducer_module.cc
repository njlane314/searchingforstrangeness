#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"

#include "Common/PandoraUtilities.h"
#include "Imaging/Image.h"
#include "Imaging/ImageCentering.h"
#include "Imaging/ImageProduction.h"
#include "Imaging/SemanticClassifier.h"
#include "Products/ImageProducts.h"

#include <TVector3.h>
#include <cetlib_except/exception.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <ios>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using image::Image;
using image::ImageProduct;
using image::ImageProperties;

class ImageProducer : public art::EDProducer {
  public:
    explicit ImageProducer(fhicl::ParameterSet const &pset);
    void produce(art::Event &e) override;

  private:
    art::InputTag fPFPproducer;
    art::InputTag fSLCproducer;
    art::InputTag fHITproducer;
    art::InputTag fWIREproducer;
    art::InputTag fMCPproducer;
    art::InputTag fBKTproducer;
    art::InputTag fSPproducer;

    bool fIsData{false};

    std::string fBadChannelFile;
    std::set<unsigned int> fBadChannels;

    std::unique_ptr<image::SemanticClassifier> fSemantic;

    int fImgW{1024};
    int fImgH{1024};
    const geo::GeometryCore *fGeo{nullptr};
    const detinfo::DetectorProperties *fDetp{nullptr};
    double fDriftStep{0.0};
    double fPitchU{0.0};
    double fPitchV{0.0};
    double fPitchW{0.0};

    void loadBadChannels(const std::string &filename);
    std::vector<art::Ptr<recob::Hit>> collectNeutrinoSliceHits(const art::Event &event) const;
    TVector3 computeImageCenter(const art::Event &event, const std::vector<art::Ptr<recob::Hit>> &neutrino_hits,
                                bool &center_defaulted) const;
};

ImageProducer::ImageProducer(fhicl::ParameterSet const &pset) {
    fPFPproducer = pset.get<art::InputTag>("PFPproducer");
    fSLCproducer = pset.get<art::InputTag>("SLCproducer");
    fHITproducer = pset.get<art::InputTag>("HITproducer");
    fWIREproducer = pset.get<art::InputTag>("WIREproducer");
    fMCPproducer = pset.get<art::InputTag>("MCPproducer");
    fBKTproducer = pset.get<art::InputTag>("BKTproducer");
    fSPproducer = pset.get<art::InputTag>("SPproducer");

    fIsData = pset.get<bool>("IsData", false);

    fBadChannelFile = pset.get<std::string>("BadChannelFile", "");
    if (!fBadChannelFile.empty())
        loadBadChannels(fBadChannelFile);

    fGeo = art::ServiceHandle<geo::Geometry>()->provider();
    fDetp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();

    auto const *det_prop_data = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->provider();
    auto const *clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->provider();

    double const tick_period = clock_data->TPCClock().TickPeriod();
    double const drift_velocity = det_prop_data->DriftVelocity();
    fDriftStep = tick_period * drift_velocity * 1.0e1;
    fPitchU = fGeo->WirePitch(geo::kU);
    fPitchV = fGeo->WirePitch(geo::kV);
    fPitchW = fGeo->WirePitch(geo::kW);

    fSemantic = std::make_unique<image::SemanticClassifier>(fMCPproducer);

    produces<std::vector<ImageProduct>>("NuSlice");
}

void ImageProducer::loadBadChannels(const std::string &filename) {
    fBadChannels.clear();
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw art::Exception(art::errors::Configuration) << "Cannot open bad channel file: " << filename;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::stringstream ss(line);
        unsigned first = 0;
        unsigned second = 0;
        ss >> first;
        if (ss >> second) {
            for (unsigned ch = first; ch <= second; ++ch)
                fBadChannels.insert(ch);
        } else {
            fBadChannels.insert(first);
        }
    }
}

std::vector<art::Ptr<recob::Hit>> ImageProducer::collectNeutrinoSliceHits(const art::Event &event) const {
    auto pfp_handle = event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::FindManyP<recob::Slice> pfp_to_slice(pfp_handle, event, fPFPproducer);

    std::optional<size_t> nu_index;
    for (size_t i = 0; i < pfp_handle->size(); ++i) {
        auto const &p = pfp_handle->at(i);
        if (!p.IsPrimary())
            continue;
        int pdg = std::abs(p.PdgCode());
        if (pdg == 12 || pdg == 14 || pdg == 16) {
            nu_index = i;
            break;
        }
    }

    std::vector<art::Ptr<recob::Hit>> out;
    auto slice_handle = event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    art::FindManyP<recob::Hit> slice_to_hits(slice_handle, event, fSLCproducer);

    if (nu_index) {
        auto slices = pfp_to_slice.at(*nu_index);
        if (!slices.empty()) {
            auto const &sl = slices.front();
            auto hits = slice_to_hits.at(sl.key());
            out.insert(out.end(), hits.begin(), hits.end());
            return out;
        }
    }

    return out;
}

TVector3 ImageProducer::computeImageCenter(const art::Event &event,
                                           const std::vector<art::Ptr<recob::Hit>> &neutrino_hits,
                                           bool &center_defaulted) const {
    std::map<art::Ptr<recob::SpacePoint>, double> sp_charge;

    auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);

    art::FindManyP<recob::SpacePoint> hit_to_sp(hit_handle, event, fSPproducer);

    for (auto const &h : neutrino_hits) {
        if (!h)
            continue;

        double q = std::max(0.f, h->Integral());
        if (!(q > 0.0))
            continue;

        auto const &sps_for_hit = hit_to_sp.at(h.key());
        if (sps_for_hit.empty())
            continue;

        double q_each = q / static_cast<double>(sps_for_hit.size());

        for (auto const &sp : sps_for_hit) {
            if (!sp)
                continue;
            sp_charge[sp] += q_each;
        }
    }

    std::vector<art::Ptr<recob::SpacePoint>> spacepoints;
    std::vector<double> weights;
    spacepoints.reserve(sp_charge.size());
    weights.reserve(sp_charge.size());

    for (auto const &kv : sp_charge) {
        spacepoints.push_back(kv.first);
        weights.push_back(kv.second);
    }

    double center_radius_max = 0.5 * std::min(fImgH * fDriftStep, fImgW * fPitchW);

    TVector3 center_world(0., 0., 0.);
    bool center_from_spacepoints = false;

    if (!spacepoints.empty()) {
        TVector3 c0(0., 0., 0.);
        double W0 = 0.0;

        const std::size_t nsp = std::min(spacepoints.size(), weights.size());

        for (std::size_t i = 0; i < nsp; ++i) {
            auto const &sp = spacepoints[i];
            double w = weights[i];
            if (!sp || !(w > 0.0))
                continue;

            auto const *xyz = sp->XYZ();
            TVector3 p(xyz[0], xyz[1], xyz[2]);

            if (!std::isfinite(p.X()) || !std::isfinite(p.Y()) || !std::isfinite(p.Z()))
                continue;

            c0 += w * p;
            W0 += w;
        }

        if (W0 > 0.0) {
            c0 *= (1.0 / W0);

            struct RadiusWeight {
                double r;
                double w;
            };

            std::vector<RadiusWeight> rw;
            rw.reserve(nsp);

            double Wtot = 0.0;

            for (std::size_t i = 0; i < nsp; ++i) {
                auto const &sp = spacepoints[i];
                double w = weights[i];
                if (!sp || !(w > 0.0))
                    continue;

                auto const *xyz = sp->XYZ();
                TVector3 p(xyz[0], xyz[1], xyz[2]);

                if (!std::isfinite(p.X()) || !std::isfinite(p.Y()) || !std::isfinite(p.Z()))
                    continue;

                double r = (p - c0).Mag();
                rw.push_back({r, w});
                Wtot += w;
            }

            if (!rw.empty() && Wtot > 0.0) {
                std::sort(rw.begin(), rw.end(), [](RadiusWeight const &a, RadiusWeight const &b) { return a.r < b.r; });

                double frac_core = 0.7;

                double target = frac_core * Wtot;
                double acc = 0.0;
                double Rc = 0.0;

                for (auto const &x : rw) {
                    acc += x.w;
                    Rc = x.r;
                    if (acc >= target)
                        break;
                }

                if (!(Rc > 0.0) || !std::isfinite(Rc))
                    Rc = center_radius_max;
                else
                    Rc = std::min(Rc, center_radius_max);

                center_world = image::trimmedCentroid3D(spacepoints, weights, c0, Rc);
                center_from_spacepoints = true;
            } else {
                center_world = c0;
                center_from_spacepoints = true;
            }
        }
    }

    center_defaulted = !center_from_spacepoints;

    if (!std::isfinite(center_world.X()) || !std::isfinite(center_world.Y()) || !std::isfinite(center_world.Z())) {
        center_world.SetXYZ(0., 0., 0.);
        center_defaulted = true;
    }

    return center_world;
}

void ImageProducer::produce(art::Event &event) {
    mf::LogDebug("ImageProducer") << "Starting image production for event " << event.id();

    auto neutrino_hits = collectNeutrinoSliceHits(event);

    if (!fBadChannels.empty()) {
        auto remove_bad_channels = [&](std::vector<art::Ptr<recob::Hit>> &hits) {
            hits.erase(std::remove_if(hits.begin(), hits.end(),
                                      [&](auto const &h) {
                                          if (h.isNull())
                                              return false;
                                          return fBadChannels.count(h->Channel()) > 0;
                                      }),
                       hits.end());
        };
        remove_bad_channels(neutrino_hits);
    }

    if (neutrino_hits.empty()) {
        mf::LogDebug("ImageProducer") << "No neutrino slice hits found; producing empty NuSlice for event "
                                      << event.id();

        auto out_slice = std::make_unique<std::vector<ImageProduct>>();
        event.put(std::move(out_slice), "NuSlice");

        mf::LogInfo("ImageProducer") << "Stored 0 NuSlice ImageProducts for event " << event.id()
                                     << " (no neutrino slice)";
        return;
    }

    bool center_defaulted = false;
    TVector3 center_world = computeImageCenter(event, neutrino_hits, center_defaulted);

    std::cout << "[ImageProducer] center_world = (" << center_world.X() << ", " << center_world.Y() << ", "
              << center_world.Z() << ") - " << (center_defaulted ? "DEFAULTED" : "COMPUTED") << std::endl;

    TVector3 cU = common::ProjectToWireView(center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_U);
    TVector3 cV = common::ProjectToWireView(center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_V);
    TVector3 cW = common::ProjectToWireView(center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_W);

    std::vector<ImageProperties> props;
    props.emplace_back(cU.Z(), cU.X(), fImgW, fImgH, fDriftStep, fPitchU, geo::kU);
    props.emplace_back(cV.Z(), cV.X(), fImgW, fImgH, fDriftStep, fPitchV, geo::kV);
    props.emplace_back(cW.Z(), cW.X(), fImgW, fImgH, fDriftStep, fPitchW, geo::kW);

    std::vector<Image<float>> det_slice;
    std::vector<Image<int>> sem_slice;

    image::PixelImageOptions opts;
    opts.producers = {fWIREproducer, fHITproducer, fMCPproducer, fBKTproducer};
    opts.semantic = fIsData ? nullptr : fSemantic.get();

    mf::LogDebug("ImageProducer") << "IsData=" << std::boolalpha << fIsData << ", semantic images "
                                  << (fIsData ? "DISABLED" : "ENABLED");

    image::ImageProduction builder(*fGeo, opts);

    builder.build(event, neutrino_hits, props, det_slice, sem_slice, fDetp);

    mf::LogDebug("ImageProducer") << "Built images: det=" << det_slice.size() << ", sem=" << sem_slice.size();

    for (std::size_t i = 0; i < props.size(); ++i) {
        mf::LogDebug("ImageProducer") << "View " << static_cast<int>(props[i].view()) << " size " << props[i].width()
                                      << "x" << props[i].height();
    }

    auto pack_plane = [](Image<float> const &det, Image<int> const &sem, ImageProperties const &p, bool include_sem) {
        ImageProduct out;
        out.view = static_cast<int>(p.view());
        out.width = static_cast<uint32_t>(p.width());
        out.height = static_cast<uint32_t>(p.height());
        out.origin_x = static_cast<float>(p.origin_x());
        out.origin_y = static_cast<float>(p.origin_y());
        out.pixel_w = static_cast<float>(p.pixel_w());
        out.pixel_h = static_cast<float>(p.pixel_h());
        out.adc = det.data();
        if (include_sem) {
            auto tmp = sem.data();
            out.semantic.assign(tmp.begin(), tmp.end());
        }
        return out;
    };

    auto out_slice = std::make_unique<std::vector<ImageProduct>>();
    out_slice->reserve(3);
    for (std::size_t i = 0; i < 3 && i < det_slice.size(); ++i) {
        out_slice->emplace_back(pack_plane(det_slice[i], sem_slice[i], props[i], !fIsData));
    }

    auto const n_slice = out_slice->size();
    event.put(std::move(out_slice), "NuSlice");

    mf::LogInfo("ImageProducer") << "Stored " << n_slice << " NuSlice ImageProducts for event " << event.id();
}

DEFINE_ART_MODULE(ImageProducer)
