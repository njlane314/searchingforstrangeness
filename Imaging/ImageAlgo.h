#ifndef IMAGEALGO_H
#define IMAGEALGO_H

#include "Imaging/Image.h"
#include "Imaging/ImageProducer.h"
#include "Imaging/InferenceEngine.h"
#include "Imaging/SemanticPixelClassifier.h"

#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <cstdint>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace analysis {

struct BinaryInferenceOutput {
    std::unordered_map<std::string, float> scores;
    bool has_segmentation{false};
    uint32_t segW{0}, segH{0};
    std::vector<uint8_t> seg_u, seg_v, seg_w;
    std::vector<float>   seg_conf_u, seg_conf_v, seg_conf_w;
};

struct ModelConfig {
    std::string name;
    std::string weights_file;
    std::string arch;
};

class ImageAlgo {
  public:
    ImageAlgo(const art::InputTag &wireProducer,
              const art::InputTag &hitProducer,
              const art::InputTag &mcpProducer,
              const art::InputTag &bktProducer,
              float adcThreshold,
              const std::string &weightsBaseDir,
              const std::string &inferenceWrapper,
              const std::string &assetsBaseDir,
              const std::vector<ModelConfig> &models,
              const std::vector<std::string> &activeModels,
              const geo::GeometryCore *geo,
              const detinfo::DetectorProperties *detp,
              const std::string &workDir);

    void produceImages(const art::Event &event,
                       const std::vector<art::Ptr<recob::Hit>> &hits,
                       const std::vector<ImageProperties> &properties,
                       bool is_data,
                       SemanticPixelClassifier *semantic_classifier,
                       const std::set<unsigned int> &badChannels,
                       std::vector<Image<float>> &detector_images,
                       std::vector<Image<int>> &semantic_images) const;

    std::unordered_map<std::string, float>
    runInference(const std::vector<Image<float>> &detector_images,
                 const std::string &absolute_scratch_dir) const;

    BinaryInferenceOutput
    runInferenceBinary(const std::vector<Image<float>> &detector_images,
                       const std::string &absolute_scratch_dir,
                       bool want_cls = true, bool want_seg = true) const;

  private:
    art::InputTag fWIREproducer;
    art::InputTag fHITproducer;
    art::InputTag fMCPproducer;
    art::InputTag fBKTproducer;
    float fADCImageThreshold;
    std::string fWeightsBaseDir;
    std::string fInferenceWrapper;
    std::string fAssetsBaseDir;
    std::vector<ModelConfig> fModels;
    std::vector<std::string> fActiveModels;
    const geo::GeometryCore *fGeo;
    const detinfo::DetectorProperties *fDetp;
    std::string fWorkDir;
};

inline bool isAbsPath(const std::string &p) { return !p.empty() && p.front() == '/'; }

inline ImageAlgo::ImageAlgo(const art::InputTag &wireProducer,
                            const art::InputTag &hitProducer,
                            const art::InputTag &mcpProducer,
                            const art::InputTag &bktProducer,
                            float adcThreshold,
                            const std::string &weightsBaseDir,
                            const std::string &inferenceWrapper,
                            const std::string &assetsBaseDir,
                            const std::vector<ModelConfig> &models,
                            const std::vector<std::string> &activeModels,
                            const geo::GeometryCore *geo,
                            const detinfo::DetectorProperties *detp,
                            const std::string &workDir)
    : fWIREproducer(wireProducer),
      fHITproducer(hitProducer),
      fMCPproducer(mcpProducer),
      fBKTproducer(bktProducer),
      fADCImageThreshold(adcThreshold),
      fWeightsBaseDir(weightsBaseDir),
      fInferenceWrapper(inferenceWrapper),
      fAssetsBaseDir(assetsBaseDir),
      fModels(models),
      fActiveModels(activeModels),
      fGeo(geo),
      fDetp(detp),
      fWorkDir(workDir) {}

inline void ImageAlgo::produceImages(
    const art::Event &event,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<ImageProperties> &properties,
    bool is_data,
    SemanticPixelClassifier *semantic_classifier,
    const std::set<unsigned int> &badChannels,
    std::vector<Image<float>> &detector_images,
    std::vector<Image<int>> &semantic_images) const {
    ImageProducer::constructPixelImages(event, hits, properties, detector_images,
                                        semantic_images, is_data, fWIREproducer,
                                        fHITproducer, fMCPproducer, fBKTproducer,
                                        fGeo, fDetp, fADCImageThreshold,
                                        semantic_classifier, badChannels);
}

inline std::unordered_map<std::string, float> ImageAlgo::runInference(
    const std::vector<Image<float>> &detector_images,
    const std::string &absolute_scratch_dir) const {
    std::unordered_map<std::string, float> scores;
    std::vector<ModelConfig> todo;
    if (fActiveModels.empty())
        todo = fModels;
    else {
        std::set<std::string> want(fActiveModels.begin(), fActiveModels.end());
        for (auto const &m : fModels)
            if (want.count(m.name))
                todo.push_back(m);
    }
    for (auto const &m : todo) {
        std::string wpath = isAbsPath(m.weights_file)
                                ? m.weights_file
                                : joinPath(fWeightsBaseDir, m.weights_file);
        float score = InferenceEngine::runInference(detector_images,
                                                    absolute_scratch_dir,
                                                    fWorkDir, m.arch, wpath,
                                                    fInferenceWrapper,
                                                    fAssetsBaseDir);
        scores[m.name] = score;
    }
    return scores;
}

inline BinaryInferenceOutput ImageAlgo::runInferenceBinary(
    const std::vector<Image<float>> &detector_images,
    const std::string &absolute_scratch_dir,
    bool want_cls, bool want_seg) const
{
    BinaryInferenceOutput out;

    std::vector<ModelConfig> todo;
    if (fActiveModels.empty())
        todo = fModels;
    else {
        std::set<std::string> want(fActiveModels.begin(), fActiveModels.end());
        for (auto const &m : fModels)
            if (want.count(m.name))
                todo.push_back(m);
    }

    bool took_seg = false;
    for (auto const& m : todo) {
        std::string wpath = isAbsPath(m.weights_file)
                                ? m.weights_file
                                : joinPath(fWeightsBaseDir, m.weights_file);
        bool ask_seg = want_seg && !took_seg;
        auto r = InferenceEngine::runInferenceDetailed(detector_images,
                                                       absolute_scratch_dir,
                                                       fWorkDir, m.arch, wpath,
                                                       fInferenceWrapper, fAssetsBaseDir,
                                                       want_cls,
                                                       ask_seg);
        if (want_cls && !r.cls.empty()) {
            out.scores[m.name] = r.cls.front();
        }
        if (ask_seg && r.segW && r.segH && !took_seg) {
            out.has_segmentation = true; out.segW = r.segW; out.segH = r.segH;
            out.seg_u = std::move(r.seg_u); out.seg_v = std::move(r.seg_v); out.seg_w = std::move(r.seg_w);
            out.seg_conf_u = std::move(r.conf_u); out.seg_conf_v = std::move(r.conf_v); out.seg_conf_w = std::move(r.conf_w);
            took_seg = true;
        }
    }
    return out;
}

}

#endif
