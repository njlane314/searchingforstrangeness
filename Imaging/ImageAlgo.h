#ifndef IMAGEALGO_H
#define IMAGEALGO_H

#include "Imaging/Image.h"
#include "Imaging/ImageProducer.h"
#include "Imaging/SemanticClassifier.h"

#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <cstdint>
#include <set>
#include <string>
#include <vector>

namespace image {

class ImageAlgo {
  public:
    ImageAlgo(const art::InputTag &wire_producer,
              const art::InputTag &hit_producer,
              const art::InputTag &mcp_producer,
              const art::InputTag &bkt_producer,
              float adc_threshold,
              const geo::GeometryCore *geo,
              const detinfo::DetectorProperties *detp);

    void produceImages(const art::Event &event,
                       const std::vector<art::Ptr<recob::Hit>> &hits,
                       const std::vector<ImageProperties> &properties,
                       bool is_data,
                       SemanticClassifier *semantic_classifier,
                       const std::set<unsigned int> &bad_channels,
                       std::vector<Image<float>> &detector_images,
                       std::vector<Image<int>> &semantic_images) const;

  private:
    art::InputTag fWIREproducer;
    art::InputTag fHITproducer;
    art::InputTag fMCPproducer;
    art::InputTag fBKTproducer;
    float fADCImageThreshold;
    const geo::GeometryCore *fGeo;
    const detinfo::DetectorProperties *fDetp;
};

inline bool isAbsPath(const std::string &p) { return !p.empty() && p.front() == '/'; }

inline ImageAlgo::ImageAlgo(const art::InputTag &wire_producer,
                            const art::InputTag &hit_producer,
                            const art::InputTag &mcp_producer,
                            const art::InputTag &bkt_producer,
                            float adc_threshold,
                            const geo::GeometryCore *geo,
                            const detinfo::DetectorProperties *detp)
    : fWIREproducer(wire_producer),
      fHITproducer(hit_producer),
      fMCPproducer(mcp_producer),
      fBKTproducer(bkt_producer),
      fADCImageThreshold(adc_threshold),
      fGeo(geo),
      fDetp(detp) {}

inline void ImageAlgo::produceImages(
    const art::Event &event,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<ImageProperties> &properties,
    bool is_data,
    SemanticClassifier *semantic_classifier,
    const std::set<unsigned int> &bad_channels,
    std::vector<Image<float>> &detector_images,
    std::vector<Image<int>> &semantic_images) const {
    ImageProducer::constructPixelImages(event, 
                                        hits, 
                                        properties, 
                                        detector_images,
                                        semantic_images, 
                                        is_data, 
                                        fWIREproducer,
                                        fHITproducer, 
                                        fMCPproducer, 
                                        fBKTproducer,
                                        fGeo, 
                                        fDetp, 
                                        fADCImageThreshold,
                                        semantic_classifier, 
                                        bad_channels);
}

}

#endif
