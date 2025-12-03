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
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"

#include "Common/PandoraUtilities.h"
#include "Imaging/Image.h"
#include "Imaging/ImageCentering.h"
#include "Imaging/ImageProduction.h"
#include "Imaging/SemanticClassifier.h"
#include "Products/ImageProducts.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <TTree.h>
#include <TVector3.h>
#include <algorithm>
#include <cmath>
#include <cetlib_except/exception.h>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <ios>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using image::Image;
using image::ImageProperties;
using image::ImageProduct;

class ImageProducer : public art::EDProducer {
public:
  explicit ImageProducer(fhicl::ParameterSet const &pset);
  void beginJob() override;
  void produce(art::Event &e) override;

private:
  art::InputTag fPFPproducer;
  art::InputTag fSLCproducer;
  art::InputTag fHITproducer;
  art::InputTag fWIREproducer;
  art::InputTag fMCPproducer;
  art::InputTag fBKTproducer;

  bool fIsData{false};

  std::string fBadChannelFile;
  std::set<unsigned int> fBadChannels;

  std::unique_ptr<sem::SemanticClassifier> fSemantic;

  int   fImgW{1024};
  int   fImgH{1024};
  float fADCThresh{4.0f};

  const geo::GeometryCore           *fGeo{nullptr};
  const detinfo::DetectorProperties *fDetp{nullptr};
  double fDriftStep{0.0};
  double fPitchU{0.0};
  double fPitchV{0.0};
  double fPitchW{0.0};

  bool        fDumpImages{false};
  std::string fDumpTreeName{"ImageDump"};

  TTree*      fImageDumpTree{nullptr};

  int   fDumpRun{0};
  int   fDumpSubrun{0};
  int   fDumpEvent{0};
  int   fDumpView{0};
  int   fDumpWidth{0};
  int   fDumpHeight{0};
  float fDumpOriginX{0.f}, fDumpOriginY{0.f};
  float fDumpPixelW{0.f},  fDumpPixelH{0.f};
  std::vector<float> fDumpADC;

  void loadBadChannels(const std::string &filename);
  std::vector<art::Ptr<recob::Hit>> collectNeutrinoSliceHits(const art::Event &event) const;
  std::vector<art::Ptr<recob::Hit>> collectEventHits(const art::Event &event) const;
};

ImageProducer::ImageProducer(fhicl::ParameterSet const &pset) {
  fPFPproducer = pset.get<art::InputTag>("PFPproducer");
  fSLCproducer = pset.get<art::InputTag>("SLCproducer");
  fHITproducer = pset.get<art::InputTag>("HITproducer");
  fWIREproducer = pset.get<art::InputTag>("WIREproducer");
  fMCPproducer = pset.get<art::InputTag>("MCPproducer");
  fBKTproducer = pset.get<art::InputTag>("BKTproducer");

  fIsData = pset.get<bool>("IsData", false);

  fBadChannelFile = pset.get<std::string>("BadChannelFile", "");
  if (!fBadChannelFile.empty())
    loadBadChannels(fBadChannelFile);

  fImgW      = pset.get<int>("ImageWidth", 1024);
  fImgH      = pset.get<int>("ImageHeight", 1024);
  fADCThresh = pset.get<float>("ADCImageThreshold", 4.0);

  fGeo  = art::ServiceHandle<geo::Geometry>()->provider();
  fDetp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();

  auto const* det_prop_data =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->provider();
  auto const* clock_data =
    art::ServiceHandle<detinfo::DetectorClocksService const>()->provider();

  double const tick_period    = clock_data->TPCClock().TickPeriod();
  double const drift_velocity = det_prop_data->DriftVelocity();
  fDriftStep = tick_period * drift_velocity * 1.0e1;
  fPitchU    = fGeo->WirePitch(geo::kU);
  fPitchV    = fGeo->WirePitch(geo::kV);
  fPitchW    = fGeo->WirePitch(geo::kW);

  fSemantic = std::make_unique<sem::SemanticClassifier>(fMCPproducer);

  fDumpImages   = pset.get<bool>("DumpImages", false);
  fDumpTreeName = pset.get<std::string>("DumpTreeName", "ImageDump");

  produces<std::vector<ImageProduct>>("NuSlice");
}

void ImageProducer::beginJob()
{
    if (!fDumpImages) return;

    auto& tfs = *art::ServiceHandle<art::TFileService>{};

    fImageDumpTree = tfs.make<TTree>(fDumpTreeName.c_str(), "Image dump");

    fImageDumpTree->Branch("run",     &fDumpRun,    "run/I");
    fImageDumpTree->Branch("subrun",  &fDumpSubrun, "subrun/I");
    fImageDumpTree->Branch("event",   &fDumpEvent,  "event/I");

    fImageDumpTree->Branch("view",    &fDumpView,   "view/I");
    fImageDumpTree->Branch("width",   &fDumpWidth,  "width/I");
    fImageDumpTree->Branch("height",  &fDumpHeight, "height/I");

    fImageDumpTree->Branch("origin_x", &fDumpOriginX, "origin_x/F");
    fImageDumpTree->Branch("origin_y", &fDumpOriginY, "origin_y/F");
    fImageDumpTree->Branch("pixel_w",  &fDumpPixelW,  "pixel_w/F");
    fImageDumpTree->Branch("pixel_h",  &fDumpPixelH,  "pixel_h/F");

    fImageDumpTree->Branch("image", &fDumpADC);
}

void ImageProducer::loadBadChannels(const std::string &filename) {
    fBadChannels.clear();
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw art::Exception(art::errors::Configuration)
            << "Cannot open bad channel file: " << filename;
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

std::vector<art::Ptr<recob::Hit>>
ImageProducer::collectNeutrinoSliceHits(const art::Event &event) const {
    auto pfp_handle =
        event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
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
    auto slice_handle =
        event.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
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

std::vector<art::Ptr<recob::Hit>>
ImageProducer::collectEventHits(const art::Event &event) const {
    auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);

    std::vector<art::Ptr<recob::Hit>> hits;
    hits.reserve(hit_handle->size());

    for (std::size_t i = 0; i < hit_handle->size(); ++i) {
        hits.emplace_back(hit_handle, i);
    }

    return hits;
}

void ImageProducer::produce(art::Event &event) {
  mf::LogDebug("ImageProducer")
    << "Starting image production for event " << event.id();

  auto neutrino_hits = collectNeutrinoSliceHits(event);
  auto event_hits    = collectEventHits(event);

  // Remove hits on bad channels, if configured.
  if (!fBadChannels.empty()) {
    auto remove_bad_channels = [&](std::vector<art::Ptr<recob::Hit>> &hits) {
      hits.erase(std::remove_if(hits.begin(), hits.end(),
                                [&](auto const &h) {
                                  if (h.isNull()) return false;
                                  return fBadChannels.count(h->Channel()) > 0;
                                }),
                 hits.end());
    };
    remove_bad_channels(neutrino_hits);
    remove_bad_channels(event_hits);
  }

  double vtx_x = std::numeric_limits<double>::quiet_NaN();
  double vtx_y = std::numeric_limits<double>::quiet_NaN();
  double vtx_z = std::numeric_limits<double>::quiet_NaN();
  {
    auto pfp_handle =
      event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::FindManyP<recob::Vertex> pfpToVtx(pfp_handle, event, fPFPproducer);
    for (std::size_t i = 0; i < pfp_handle->size(); ++i) {
      auto const &p = pfp_handle->at(i);
      if (!p.IsPrimary())
        continue;
      auto const &v = pfpToVtx.at(i);
      if (!v.empty()) {
        auto const &pos = v.front()->position();
        vtx_x = pos.X();
        vtx_y = pos.Y();
        vtx_z = pos.Z();
        break;
      }
    }
  }

  TVector3 vtx_world(vtx_x, vtx_y, vtx_z);

  TVector3 center_world = image::robustGlobalCentroid(
    event, neutrino_hits, static_cast<double>(fADCThresh), 2.0);

  if (!std::isfinite(center_world.X()) ||
      !std::isfinite(center_world.Y()) ||
      !std::isfinite(center_world.Z())) {
    center_world = vtx_world;
  }

  // Project centre to each view
  TVector3 cU = common::ProjectToWireView(
    center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_U);
  TVector3 cV = common::ProjectToWireView(
    center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_V);
  TVector3 cW = common::ProjectToWireView(
    center_world.X(), center_world.Y(), center_world.Z(), common::TPC_VIEW_W);

  std::vector<ImageProperties> props;
  props.emplace_back(cU.Z(), cU.X(),
                     fImgW, fImgH, fDriftStep,
                     fPitchU, geo::kU);
  props.emplace_back(cV.Z(), cV.X(),
                     fImgW, fImgH, fDriftStep,
                     fPitchV, geo::kV);
  props.emplace_back(cW.Z(), cW.X(),
                     fImgW, fImgH, fDriftStep,
                     fPitchW, geo::kW);

  std::vector<Image<float>> det_slice;
  std::vector<Image<int>>   sem_slice;

  image::PixelImageOptions opts;
  opts.producers = {fWIREproducer, fHITproducer, fMCPproducer, fBKTproducer};
  opts.semantic  = fIsData ? nullptr : fSemantic.get();

  mf::LogDebug("ImageProducer")
    << "IsData=" << std::boolalpha << fIsData
    << ", semantic images " << (fIsData ? "DISABLED" : "ENABLED");

  image::ImageProduction builder(*fGeo, opts);

  builder.build(
    event,
    event_hits,
    props,
    det_slice,
    sem_slice,
    fDetp);

  mf::LogDebug("ImageProducer")
    << "Built images: det=" << det_slice.size()
    << ", sem=" << sem_slice.size();

  for (std::size_t i = 0; i < props.size(); ++i) {
    mf::LogDebug("ImageProducer")
      << "View " << static_cast<int>(props[i].view())
      << " size " << props[i].width() << "x" << props[i].height();
  }

  // Optional dump to ROOT TTree
  if (fDumpImages && fImageDumpTree) {
    fDumpRun    = event.id().run();
    fDumpSubrun = event.id().subRun();
    fDumpEvent  = event.id().event();

    std::vector<Image<float>> const* images = &det_slice;

    for (std::size_t i = 0; i < images->size() && i < props.size(); ++i) {
      auto const& img  = images->at(i);
      auto const& prop = props[i];

      fDumpView    = static_cast<int>(prop.view());
      fDumpWidth   = static_cast<int>(prop.width());
      fDumpHeight  = static_cast<int>(prop.height());
      fDumpOriginX = static_cast<float>(prop.origin_x());
      fDumpOriginY = static_cast<float>(prop.origin_y());
      fDumpPixelW  = static_cast<float>(prop.pixel_w());
      fDumpPixelH  = static_cast<float>(prop.pixel_h());

      fDumpADC = img.data();

      fImageDumpTree->Fill();
    }
  }

  // Pack into ImageProduct and store in event
  auto pack_plane = [](Image<float> const &det, Image<int> const &sem,
                       ImageProperties const &p, bool include_sem) {
    ImageProduct out;
    out.view    = static_cast<int>(p.view());
    out.width   = static_cast<uint32_t>(p.width());
    out.height  = static_cast<uint32_t>(p.height());
    out.origin_x = static_cast<float>(p.origin_x());
    out.origin_y = static_cast<float>(p.origin_y());
    out.pixel_w  = static_cast<float>(p.pixel_w());
    out.pixel_h  = static_cast<float>(p.pixel_h());
    out.adc      = det.data();
    if (include_sem) {
      auto tmp = sem.data();
      out.semantic.assign(tmp.begin(), tmp.end());
    }
    return out;
  };

  auto out_slice = std::make_unique<std::vector<ImageProduct>>();
  out_slice->reserve(3);
  for (std::size_t i = 0; i < 3 && i < det_slice.size(); ++i) {
    out_slice->emplace_back(
      pack_plane(det_slice[i], sem_slice[i], props[i], !fIsData));
  }

  auto const n_slice = out_slice->size();
  event.put(std::move(out_slice), "NuSlice");

  mf::LogInfo("ImageProducer")
    << "Stored " << n_slice
    << " NuSlice ImageProducts for event " << event.id();
}

DEFINE_ART_MODULE(ImageProducer)
