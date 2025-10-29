#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "Imaging/Image.h"
#include "Imaging/InferenceProduction.h"
#include "Products/ImageProducts.h"
#include "Products/InferencePerf.h"
#include "Products/RandomFeatures.h"
#include "Products/SegmentationProducts.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>
#include <system_error>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

std::string join_path(std::string base, const std::string &leaf) {
  if (leaf.empty()) return base;
  if (base.empty()) return leaf;
  if (base.back() != '/') base.push_back('/');
  return base + leaf;
}

std::string resolve_base_dir(const std::string &configured, const char *env_var) {
  if (env_var && *env_var) {
    if (const char *env = std::getenv(env_var)) {
      if (*env) return env;
    }
  }
  return configured;
}

std::string make_absolute_path(const std::string &path) {
  if (path.empty()) return std::filesystem::current_path().string();
  std::error_code ec;
  auto abs = std::filesystem::absolute(path, ec);
  return ec ? path : abs.string();
}

std::string resolve_with_base(const std::string &base, const std::string &path) {
  if (path.empty()) return path;
  if (!path.empty() && (path.front() == '/' || (path.size() > 1 && path[1] == ':')))
    return path;
  return join_path(base, path);
}

image::Image<float> plane_to_image(const image::PlaneImage &plane) {
  const double centre_x = plane.origin_x + 0.5 * plane.pixel_w * plane.width;
  const double centre_y = plane.origin_y + 0.5 * plane.pixel_h * plane.height;
  image::ImageProperties props(centre_x, centre_y,
                               static_cast<size_t>(plane.width), static_cast<size_t>(plane.height),
                               plane.pixel_h, plane.pixel_w,
                               static_cast<geo::View_t>(plane.view));
  image::Image<float> img(props);
  const auto expected = static_cast<size_t>(plane.width) * static_cast<size_t>(plane.height);
  if (plane.adc.size() != expected) {
    throw cet::exception("InferenceProduction")
      << "Plane ADC vector has size " << plane.adc.size()
      << " but expected " << expected;
  }
  for (size_t row = 0; row < plane.height; ++row) {
    const size_t offset = row * plane.width;
    for (size_t col = 0; col < plane.width; ++col) {
      img.set(row, col, plane.adc[offset + col], false);
    }
  }
  return img;
}

std::vector<std::string> parse_active_models(const fhicl::ParameterSet &p) {
  std::vector<std::string> active = p.get<std::vector<std::string>>("ActiveModels", {});
  if (active.empty()) {
    if (const char *env = std::getenv("ACTIVE_MODELS")) {
      std::string from_env(env);
      std::string current;
      for (char c : from_env) {
        if (c == ',' || c == ';' || std::isspace(static_cast<unsigned char>(c))) {
          if (!current.empty()) {
            active.push_back(current);
            current.clear();
          }
        } else {
          current.push_back(c);
        }
      }
      if (!current.empty()) active.push_back(current);
    }
  }
  return active;
}

} // namespace

namespace image {

class InferenceProduction : public art::EDProducer {
public:
  explicit InferenceProduction(fhicl::ParameterSet const &p);
  void produce(art::Event &event) override;

private:
  struct ModelConfig {
    std::string name;
    std::string arch;
    std::string weights;
    std::vector<std::string> classNames;
    bool wantClassification{true};
    bool wantSegmentation{false};
    bool storeFeatures{false};
  };

  art::InputTag planesTag_;
  std::string assetsBaseDir_;
  std::string weightsBaseDir_;
  std::string inferenceWrapper_;
  std::string scratchDir_;
  bool wantClassificationDefault_{true};
  bool wantSegmentationDefault_{false};
  bool storeFeaturesDefault_{false};
  std::vector<ModelConfig> models_;

  std::string scratch_dir() const;
};

InferenceProduction::InferenceProduction(fhicl::ParameterSet const &p) {
  planesTag_ = p.get<art::InputTag>("PlanesTag");

  assetsBaseDir_ = resolve_base_dir(p.get<std::string>("AssetsBaseDir", ""), "ASSETS_BASE_DIR");
  weightsBaseDir_ = resolve_base_dir(p.get<std::string>("WeightsBaseDir", assetsBaseDir_), "WEIGHTS_BASE_DIR");

  inferenceWrapper_ = p.get<std::string>("InferenceWrapper", join_path(assetsBaseDir_, "scripts/inference_wrapper.sh"));
  if (const char *env = std::getenv("INFERENCE_WRAPPER")) {
    if (*env) inferenceWrapper_ = env;
  }

  scratchDir_ = p.get<std::string>("ScratchDir", "");
  wantClassificationDefault_ = p.get<bool>("WantClassification", true);
  wantSegmentationDefault_ = p.get<bool>("WantSegmentation", false);
  storeFeaturesDefault_ = p.get<bool>("StoreFeatures", false);

  auto activeModels = parse_active_models(p);
  std::unordered_set<std::string> activeSet(activeModels.begin(), activeModels.end());
  bool filterActive = !activeSet.empty();

  auto model_psets = p.get<std::vector<fhicl::ParameterSet>>("Models", {});
  for (auto const &mp : model_psets) {
    ModelConfig cfg;
    cfg.name = mp.get<std::string>("Name");
    if (filterActive && activeSet.count(cfg.name) == 0) continue;
    cfg.arch = mp.get<std::string>("Architecture", "");
    cfg.weights = resolve_with_base(weightsBaseDir_, mp.get<std::string>("Weights", ""));
    cfg.classNames = mp.get<std::vector<std::string>>("ClassNames", {});
    cfg.wantClassification = mp.get<bool>("WantClassification", wantClassificationDefault_);
    cfg.wantSegmentation = mp.get<bool>("WantSegmentation", wantSegmentationDefault_);
    cfg.storeFeatures = mp.get<bool>("StoreFeatures", storeFeaturesDefault_);
    models_.push_back(std::move(cfg));
  }

  if (models_.empty()) {
    mf::LogWarning("InferenceProduction")
        << "No inference models configured; module will produce empty outputs.";
  }

  produces<std::vector<image::InferenceScores>>();
  produces<std::vector<image::PlaneSegmentation>>("seg");
  produces<image::InferencePerfProduct>("perf");
  produces<std::vector<image::RandomFeatures>>("features");
}

std::string InferenceProduction::scratch_dir() const {
  if (!scratchDir_.empty()) return scratchDir_;
  if (const char *env = std::getenv("_CONDOR_SCRATCH_DIR")) return env;
  if (const char *tmp = std::getenv("TMPDIR")) return tmp;
  return ".";
}

void InferenceProduction::produce(art::Event &event) {
  auto planes = event.getValidHandle<std::vector<image::PlaneImage>>(planesTag_);
  if (planes->size() < 3) {
    throw cet::exception("InferenceProduction")
        << "Need at least three planes (U, V, W); got " << planes->size();
  }

  auto const *U = &planes->at(0);
  auto const *V = &planes->at(1);
  auto const *W = &planes->at(2);
  for (auto const &plane : *planes) {
    if      (plane.view == static_cast<int>(geo::kU)) U = &plane;
    else if (plane.view == static_cast<int>(geo::kV)) V = &plane;
    else if (plane.view == static_cast<int>(geo::kW)) W = &plane;
  }

  const uint32_t width = U->width;
  const uint32_t height = U->height;
  auto check_dims = [&](auto const *plane) {
    if (!plane) return;
    if (plane->width != width || plane->height != height) {
      throw cet::exception("InferenceProduction")
          << "Plane size mismatch. Expected " << width << "x" << height
          << ", got " << plane->width << "x" << plane->height;
    }
  };
  check_dims(V);
  check_dims(W);

  std::vector<image::Image<float>> detector_images;
  detector_images.reserve(3);
  detector_images.emplace_back(plane_to_image(*U));
  detector_images.emplace_back(plane_to_image(*V));
  detector_images.emplace_back(plane_to_image(*W));

  const std::string scratch = scratch_dir();
  const std::string absScratch = make_absolute_path(scratch);

  auto scoresOut = std::make_unique<std::vector<image::InferenceScores>>();
  auto segOut = std::make_unique<std::vector<image::PlaneSegmentation>>();
  auto perfOut = std::make_unique<image::InferencePerfProduct>();
  auto featuresOut = std::make_unique<std::vector<image::RandomFeatures>>();

  for (auto const &model : models_) {
    auto result = InferenceEngine::runInferenceDetailed(detector_images,
                                                        absScratch,
                                                        scratch,
                                                        model.arch,
                                                        model.weights,
                                                        inferenceWrapper_,
                                                        assetsBaseDir_,
                                                        model.wantClassification,
                                                        model.wantSegmentation);

    image::ModelPerf perf{};
    perf.model = model.name;
    perf.t_write_req_ms  = static_cast<float>(result.perf.t_write_req_ms);
    perf.t_exec_total_ms = static_cast<float>(result.perf.t_exec_total_ms);
    perf.t_read_resp_ms  = static_cast<float>(result.perf.t_read_resp_ms);
    perf.t_child_total_ms = static_cast<float>(result.perf.t_child_total_ms);
    perf.t_child_setup_ms = static_cast<float>(result.perf.t_child_setup_ms);
    perf.t_child_infer_ms = static_cast<float>(result.perf.t_child_infer_ms);
    perf.t_child_post_ms  = static_cast<float>(result.perf.t_child_post_ms);
    perf.child_max_rss_mb = static_cast<float>(result.perf.child_max_rss_mb);
    perf.child_cuda_mem_mb = static_cast<float>(result.perf.child_cuda_mem_mb);
    perfOut->per_model.push_back(std::move(perf));

    if (model.wantClassification) {
      image::InferenceScores scores{};
      scores.model = model.name;
      scores.scores = result.cls;
      scores.names = model.classNames;
      if (scores.names.size() < scores.scores.size()) {
        for (size_t i = scores.names.size(); i < scores.scores.size(); ++i) {
          scores.names.emplace_back("score_" + std::to_string(i));
        }
      }
      if (!scores.scores.empty() || !scores.names.empty()) {
        scoresOut->push_back(std::move(scores));
      }
    }

    if (model.wantSegmentation && result.segW && result.segH) {
      auto append_seg = [&](const image::PlaneImage *plane,
                            const std::vector<uint8_t> &labels,
                            const std::vector<float> &confidence) {
        if (!plane || labels.empty()) return;
        image::PlaneSegmentation seg{};
        seg.model = model.name;
        seg.view = plane->view;
        seg.width = result.segW;
        seg.height = result.segH;
        seg.nclasses = model.classNames.empty()
                         ? 0u
                         : static_cast<uint32_t>(model.classNames.size());
        seg.labels = labels;
        if (!confidence.empty()) seg.confidence = confidence;
        segOut->push_back(std::move(seg));
      };
      append_seg(U, result.seg_u, result.conf_u);
      append_seg(V, result.seg_v, result.conf_v);
      append_seg(W, result.seg_w, result.conf_w);
    }

    if (model.storeFeatures && !result.features.empty()) {
      image::RandomFeatures features{};
      features.model = model.name;
      features.seed = result.feature_seed;
      features.values = result.features;
      featuresOut->push_back(std::move(features));
    }
  }

  event.put(std::move(scoresOut));
  event.put(std::move(segOut), "seg");
  event.put(std::move(perfOut), "perf");
  event.put(std::move(featuresOut), "features");
}

} // namespace image

DEFINE_ART_MODULE(image::InferenceProduction)
