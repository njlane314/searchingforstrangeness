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

#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace analysis {

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

} // namespace analysis

#endif // IMAGEALGO_H
