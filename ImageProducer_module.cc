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

#include "Support/PandoraUtilities.h"
#include "ImagePipeline/Image.h"
#include "ImagePipeline/ImageCentering.h"
#include "ImagePipeline/ImageProduction.h"
#include "ImagePipeline/SemanticClassifier.h"
#include "Products/SparsePlaneImage.h"

#include <TVector3.h>
#include <cetlib_except/exception.h>

#include <algorithm>
#include <array>
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
using image::ImageProperties;
using image::SparsePlaneImage;

namespace {

common::PandoraView pandoraViewFor(geo::View_t view) {
    switch (view) {
    case geo::kU:
        return common::TPC_VIEW_U;
    case geo::kV:
        return common::TPC_VIEW_V;
    case geo::kW:
    case geo::kY:
        return common::TPC_VIEW_W;
    default:
        throw cet::exception("ImageProducer") << "Unsupported geo::View_t: " << static_cast<int>(view);
    }
}

ImageProperties makeImagePropertiesFromBounds(double min_wire_coord, double max_wire_coord, double min_drift,
                                              double max_drift, double pixel_w, double pixel_h, geo::View_t view) {
    if (!(std::isfinite(min_wire_coord) && std::isfinite(max_wire_coord) && std::isfinite(min_drift) &&
          std::isfinite(max_drift))) {
        throw cet::exception("ImageProducer") << "Cannot build image properties from non-finite bounds.";
    }
    if (!(pixel_w > 0.0) || !(pixel_h > 0.0)) {
        throw cet::exception("ImageProducer") << "Cannot build image properties with non-positive pixel size.";
    }
    if (max_wire_coord < min_wire_coord || max_drift < min_drift) {
        throw cet::exception("ImageProducer") << "Cannot build image properties from inverted bounds.";
    }

    auto const width = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil((max_wire_coord - min_wire_coord) / pixel_w)) + 1);
    auto const height = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil((max_drift - min_drift) / pixel_h)) + 1);

    double const center_wire_coord = min_wire_coord + 0.5 * static_cast<double>(width) * pixel_w;
    double const center_drift = min_drift + 0.5 * static_cast<double>(height) * pixel_h;

    return ImageProperties(center_wire_coord, center_drift, width, height, pixel_h, pixel_w, view);
}

} // namespace

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
    art::InputTag fVTXproducer;

    bool fIsData{false};

    std::string fBadChannelFile;
    std::set<unsigned int> fBadChannels;

    std::unique_ptr<image::SemanticClassifier> fSemantic;

    int fImgW{512};
    int fImgH{512};
    const geo::GeometryCore *fGeo{nullptr};
    const detinfo::DetectorProperties *fDetp{nullptr};
    double fDriftStep{0.0};
    double fPitchU{0.0};
    double fPitchV{0.0};
    double fPitchW{0.0};
    std::vector<ImageProperties> fFullWindowProps;

    void loadBadChannels(const std::string &filename);
    std::vector<art::Ptr<recob::Hit>> collectEventHits(const art::Event &event) const;
    std::vector<art::Ptr<recob::Hit>> collectNeutrinoSliceHits(const art::Event &event) const;
    std::optional<TVector3> findNeutrinoVertex(const art::Event &event) const;
    TVector3 computeImageCenter(const art::Event &event, const std::vector<art::Ptr<recob::Hit>> &neutrino_hits,
                                std::optional<TVector3> const &vertex, bool &center_defaulted) const;
    std::vector<ImageProperties> computeFullWindowProperties() const;
};

ImageProducer::ImageProducer(fhicl::ParameterSet const &pset) {
    fPFPproducer = pset.get<art::InputTag>("PFPproducer");
    fSLCproducer = pset.get<art::InputTag>("SLCproducer");
    fHITproducer = pset.get<art::InputTag>("HITproducer");
    fWIREproducer = pset.get<art::InputTag>("WIREproducer");
    fMCPproducer = pset.get<art::InputTag>("MCPproducer");
    fBKTproducer = pset.get<art::InputTag>("BKTproducer");
    fSPproducer = pset.get<art::InputTag>("SPproducer");
    fVTXproducer = pset.get<art::InputTag>("VTXproducer", fPFPproducer);

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
    fFullWindowProps = computeFullWindowProperties();

    fSemantic = std::make_unique<image::SemanticClassifier>(fMCPproducer);

    produces<std::vector<SparsePlaneImage>>("CroppedWindow");
    produces<std::vector<SparsePlaneImage>>("FullWindow");
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

std::vector<art::Ptr<recob::Hit>> ImageProducer::collectEventHits(const art::Event &event) const {
    auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);

    std::vector<art::Ptr<recob::Hit>> hits;
    hits.reserve(hit_handle->size());

    for (std::size_t i = 0; i < hit_handle->size(); ++i) {
        hits.emplace_back(hit_handle, i);
    }

    return hits;
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

std::optional<TVector3> ImageProducer::findNeutrinoVertex(const art::Event &event) const {
    auto pfp_handle = event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);

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
    if (!nu_index)
        return std::nullopt;

    try {
        art::FindManyP<recob::Vertex> pfp_to_vtx(pfp_handle, event, fVTXproducer);
        auto const vtxs = pfp_to_vtx.at(*nu_index);
        if (vtxs.empty() || !vtxs.front())
            return std::nullopt;

        double xyz[3] = {0., 0., 0.};
        vtxs.front()->XYZ(xyz);
        TVector3 v(xyz[0], xyz[1], xyz[2]);

        if (!std::isfinite(v.X()) || !std::isfinite(v.Y()) || !std::isfinite(v.Z()))
            return std::nullopt;

        return v;
    } catch (cet::exception const &ex) {
        mf::LogDebug("ImageProducer") << "PFParticle->Vertex association unavailable for VTXproducer='"
                                      << fVTXproducer.label() << "': " << ex;
    }

    return std::nullopt;
}

TVector3 ImageProducer::computeImageCenter(const art::Event &event,
                                           const std::vector<art::Ptr<recob::Hit>> &neutrino_hits,
                                           std::optional<TVector3> const &vertex,
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

    TVector3 center_world(0., 0., 0.);
    TVector3 center_reference(0., 0., 0.);
    if (vertex) {
        center_reference = *vertex;
    } else {
        center_defaulted = true;
    }

    double Rc = 0.5 * std::min(fImgH * fDriftStep, fImgW * fPitchW);
    if (Rc > 0.0 && std::isfinite(Rc) && !spacepoints.empty()) {
        center_world = image::trimmedCentroid(spacepoints, weights, center_reference, Rc);
    } else {
        center_world = center_reference;
    }

    if (!std::isfinite(center_world.X()) || !std::isfinite(center_world.Y()) || !std::isfinite(center_world.Z())) {
        center_world.SetXYZ(0., 0., 0.);
        center_defaulted = true;
    }

    return center_world;
}

std::vector<ImageProperties> ImageProducer::computeFullWindowProperties() const {
    auto const &tpc = fGeo->TPC();
    auto const bounds = tpc.ActiveBoundingBox();
    double const min_drift = bounds.MinX();
    double const max_drift = bounds.MaxX();

    auto make_view_properties = [&](unsigned int plane_idx, geo::View_t view, double pitch) {
        auto const &plane = tpc.Plane(plane_idx);

        double min_wire_coord = std::numeric_limits<double>::max();
        double max_wire_coord = std::numeric_limits<double>::lowest();
        auto const pandora_view = pandoraViewFor(view);

        for (std::size_t wire_idx = 0; wire_idx < plane.Nwires(); ++wire_idx) {
            TVector3 const center = plane.Wire(static_cast<unsigned int>(wire_idx)).GetCenter();
            TVector3 const proj = common::ProjectToWireView(center.X(), center.Y(), center.Z(), pandora_view);
            min_wire_coord = std::min(min_wire_coord, static_cast<double>(proj.Z()));
            max_wire_coord = std::max(max_wire_coord, static_cast<double>(proj.Z()));
        }

        return makeImagePropertiesFromBounds(min_wire_coord, max_wire_coord, min_drift, max_drift, pitch, fDriftStep,
                                             view);
    };

    std::vector<ImageProperties> props;
    props.reserve(3);
    props.emplace_back(make_view_properties(0U, geo::kU, fPitchU));
    props.emplace_back(make_view_properties(1U, geo::kV, fPitchV));
    props.emplace_back(make_view_properties(2U, geo::kW, fPitchW));
    return props;
}

void ImageProducer::produce(art::Event &event) {
    mf::LogDebug("ImageProducer") << "Starting image production for event " << event.id();

    auto neutrino_hits = collectNeutrinoSliceHits(event);
    auto event_hits = collectEventHits(event);

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
        remove_bad_channels(event_hits);
    }

    auto out_slice = std::make_unique<std::vector<SparsePlaneImage>>();
    auto out_full = std::make_unique<std::vector<SparsePlaneImage>>();

    image::PixelImageOptions opts;
    opts.producers = {fWIREproducer, fHITproducer, fMCPproducer, fBKTproducer};
    opts.semantic = fIsData ? nullptr : fSemantic.get();

    mf::LogDebug("ImageProducer") << "IsData=" << std::boolalpha << fIsData << ", semantic images "
                                  << (fIsData ? "DISABLED" : "ENABLED");

    image::ImageProduction builder(*fGeo, opts);

    std::array<uint32_t, 3> hit_counts{0, 0, 0};
    for (auto const &hit : event_hits) {
        switch (hit->View()) {
        case geo::kU:
            ++hit_counts[0];
            break;
        case geo::kV:
            ++hit_counts[1];
            break;
        case geo::kW:
            ++hit_counts[2];
            break;
        default:
            break;
        }
    }

    auto vertex = findNeutrinoVertex(event);

    auto vertex_pixel = [&](ImageProperties const &p) -> std::pair<int32_t, int32_t> {
        if (!vertex)
            return {-1, -1};
        auto const pandora_view = pandoraViewFor(p.view());

        auto proj = common::ProjectToWireView(vertex->X(), vertex->Y(), vertex->Z(), pandora_view);
        auto col = p.col(proj.Z());
        auto row = p.row(proj.X());
        if (!col || !row)
            return {-1, -1};
        return {static_cast<int32_t>(*row), static_cast<int32_t>(*col)};
    };

    auto pack_plane = [&](Image<float> const &det, Image<int> const &sem, Image<uint8_t> const &slice,
                          ImageProperties const &p, bool include_sem, uint32_t hit_count) {
        SparsePlaneImage out;
        out.view = static_cast<int>(p.view());
        out.width = static_cast<uint32_t>(p.width());
        out.height = static_cast<uint32_t>(p.height());
        out.hit_count = hit_count;
        out.origin_x = static_cast<float>(p.origin_x());
        out.origin_y = static_cast<float>(p.origin_y());
        out.pixel_w = static_cast<float>(p.pixel_w());
        out.pixel_h = static_cast<float>(p.pixel_h());
        auto [vertex_row, vertex_col] = vertex_pixel(p);
        out.vertex_row = vertex_row;
        out.vertex_col = vertex_col;
        out.feature_dim = 2;
        auto const &det_pixels = det.data();
        auto const &sem_pixels = sem.data();
        auto const &slice_pixels = slice.data();
        auto const width = p.width();

        auto const n_active = static_cast<std::size_t>(
            std::count_if(det_pixels.begin(), det_pixels.end(),
                          [](float a) { return a > 0.f; }));

        out.coords.reserve(n_active * 2);
        out.features.reserve(n_active * out.feature_dim);
        if (include_sem)
            out.semantic.reserve(n_active);

        for (std::size_t idx = 0; idx < det_pixels.size(); ++idx) {
            float const a = det_pixels[idx];
            if (a <= 0.f)
                continue;

            auto const row = static_cast<int32_t>(idx / width);
            auto const col = static_cast<int32_t>(idx % width);

            out.coords.push_back(row);
            out.coords.push_back(col);
            out.features.push_back(a);
            out.features.push_back(static_cast<float>(slice_pixels[idx]));
            if (include_sem)
                out.semantic.push_back(static_cast<uint8_t>(sem_pixels[idx]));
        }
        return out;
    };

    auto append_planes = [&](std::vector<SparsePlaneImage> &out, std::vector<Image<float>> const &det_images,
                             std::vector<Image<int>> const &sem_images, std::vector<Image<uint8_t>> const &slice_images,
                             std::vector<ImageProperties> const &properties) {
        out.reserve(3);
        for (std::size_t i = 0; i < 3 && i < det_images.size(); ++i) {
            out.emplace_back(pack_plane(det_images[i], sem_images[i], slice_images[i], properties[i], !fIsData,
                                        hit_counts[i]));
        }
    };

    if (neutrino_hits.empty()) {
        mf::LogDebug("ImageProducer") << "No neutrino slice hits found; CroppedWindow will be empty for event "
                                      << event.id();
    } else {
        bool center_defaulted = false;
        TVector3 center_world = computeImageCenter(event, neutrino_hits, vertex, center_defaulted);

        std::cout << "[ImageProducer] center_world = (" << center_world.X() << ", " << center_world.Y() << ", "
                  << center_world.Z() << ") - " << (center_defaulted ? "NO VERTEX (FALLBACK)" : "VERTEX-SEEDED")
                  << std::endl;

        TVector3 cU =
            common::ProjectToWireView(center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_U);
        TVector3 cV =
            common::ProjectToWireView(center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_V);
        TVector3 cW =
            common::ProjectToWireView(center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_W);

        std::vector<ImageProperties> slice_props;
        slice_props.emplace_back(cU.Z(), cU.X(), fImgW, fImgH, fDriftStep, fPitchU, geo::kU);
        slice_props.emplace_back(cV.Z(), cV.X(), fImgW, fImgH, fDriftStep, fPitchV, geo::kV);
        slice_props.emplace_back(cW.Z(), cW.X(), fImgW, fImgH, fDriftStep, fPitchW, geo::kW);

        std::vector<Image<float>> det_slice;
        std::vector<Image<int>> sem_slice;
        std::vector<Image<uint8_t>> slice_mask;
        builder.build(event, event_hits, neutrino_hits, slice_props, det_slice, sem_slice, slice_mask, fDetp);

        mf::LogDebug("ImageProducer") << "Built CroppedWindow images: det=" << det_slice.size()
                                      << ", sem=" << sem_slice.size();

        for (std::size_t i = 0; i < slice_props.size(); ++i) {
            mf::LogDebug("ImageProducer") << "CroppedWindow view " << static_cast<int>(slice_props[i].view()) << " size "
                                          << slice_props[i].width() << "x" << slice_props[i].height();
        }

        append_planes(*out_slice, det_slice, sem_slice, slice_mask, slice_props);
    }

    std::vector<Image<float>> det_full;
    std::vector<Image<int>> sem_full;
    std::vector<Image<uint8_t>> full_slice_mask;
    builder.build(event, event_hits, neutrino_hits, fFullWindowProps, det_full, sem_full, full_slice_mask, fDetp);

    mf::LogDebug("ImageProducer") << "Built FullWindow images: det=" << det_full.size()
                                  << ", sem=" << sem_full.size();

    for (std::size_t i = 0; i < fFullWindowProps.size(); ++i) {
        mf::LogDebug("ImageProducer") << "FullWindow view " << static_cast<int>(fFullWindowProps[i].view())
                                      << " size " << fFullWindowProps[i].width() << "x"
                                      << fFullWindowProps[i].height();
    }

    append_planes(*out_full, det_full, sem_full, full_slice_mask, fFullWindowProps);

    auto const n_slice = out_slice->size();
    auto const n_full = out_full->size();
    event.put(std::move(out_slice), "CroppedWindow");
    event.put(std::move(out_full), "FullWindow");

    mf::LogInfo("ImageProducer") << "Stored " << n_slice << " CroppedWindow and " << n_full
                                 << " FullWindow SparsePlaneImages for event " << event.id();
}

DEFINE_ART_MODULE(ImageProducer)
