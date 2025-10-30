#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "cetlib_except/exception.h"

#include "Imaging/InferenceProduction.h"
#include "Products/ImageProducts.h"
#include "Products/InferencePerf.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

std::string join_path(std::string base, const std::string &leaf) {
    if (leaf.empty())
        return base;
    if (base.empty())
        return leaf;
    if (base.back() != '/')
        base.push_back('/');
    return base + leaf;
}

bool is_relative(const std::string &path) {
    if (path.empty())
        return false;
    if (path.front() == '/')
        return false;
    if (path.size() > 1 && path[1] == ':')
        return false; // Windows drive
    if (path.find("://") != std::string::npos)
        return false;
    if (path.rfind("./", 0) == 0 || path.rfind("../", 0) == 0)
        return true;
    return true;
}

std::string resolve_under(const std::string &base, const std::string &value) {
    if (value.empty())
        return value;
    if (base.empty())
        return value;
    if (!is_relative(value))
        return value;
    return join_path(base, value);
}

const image::PlaneImage *find_view(std::vector<image::PlaneImage> const &planes,
                                   geo::View_t view,
                                   const image::PlaneImage *fallback) {
    for (auto const &plane : planes) {
        if (plane.view == static_cast<int>(view))
            return &plane;
    }
    return fallback;
}

} // namespace

namespace image {

class InferenceProducerModule : public art::EDProducer {
  public:
    struct ModelConfig {
        std::string label;
        std::string arch;
        std::string weights;
        std::string wrapper;
        std::string assets;
        bool want_cls{true};
        bool want_seg{false};
    };

    explicit InferenceProducerModule(fhicl::ParameterSet const &p);

    void produce(art::Event &e) override;

  private:
    static ModelConfig makeConfig(fhicl::ParameterSet const &p);

    art::InputTag planesTag_;
    std::string scratchDir_;
    std::string assetsBaseDir_;
    std::string defaultWrapper_;
    std::vector<ModelConfig> models_;
};

InferenceProducerModule::ModelConfig
InferenceProducerModule::makeConfig(fhicl::ParameterSet const &p) {
    ModelConfig cfg;
    cfg.label = p.get<std::string>("Label", "");
    cfg.arch = p.get<std::string>("Arch");
    cfg.weights = p.get<std::string>("Weights", "");
    cfg.wrapper = p.get<std::string>("Wrapper", "");
    cfg.assets = p.get<std::string>("Assets", "");
    cfg.want_cls = p.get<bool>("WantCls", true);
    cfg.want_seg = p.get<bool>("WantSeg", false);
    return cfg;
}

InferenceProducerModule::InferenceProducerModule(
    fhicl::ParameterSet const &p)
    : planesTag_{p.get<std::string>("PlanesTag")},
      scratchDir_{p.get<std::string>("ScratchDir", "")},
      assetsBaseDir_{p.get<std::string>("AssetsBaseDir", "")},
      defaultWrapper_{resolve_under(
          assetsBaseDir_, p.get<std::string>("DefaultWrapper",
                                             "scripts/inference_wrapper.sh"))} {
    produces<image::InferencePerfProduct>();

    auto modelSets = p.get<std::vector<fhicl::ParameterSet>>("Models", {});
    models_.reserve(modelSets.size());
    for (auto const &modelPs : modelSets) {
        models_.push_back(makeConfig(modelPs));
    }
}

void InferenceProducerModule::produce(art::Event &e) {
    auto planes = e.getValidHandle<std::vector<image::PlaneImage>>(planesTag_);
    if (planes->size() < 3) {
        throw cet::exception("InferenceProduction")
            << "Need at least three planes (U,V,W), got " << planes->size();
    }

    auto const *first = &planes->front();
    auto const *U = find_view(*planes, geo::kU, first);
    auto const *V = find_view(*planes, geo::kV, first);
    auto const *W = find_view(*planes, geo::kW, first);
    if (U == nullptr || V == nullptr || W == nullptr) {
        throw cet::exception("InferenceProduction")
            << "Unable to identify U/V/W planes";
    }

    std::vector<image::PlaneImage> detector_images;
    detector_images.reserve(3);
    detector_images.emplace_back(*U);
    detector_images.emplace_back(*V);
    detector_images.emplace_back(*W);

    std::string scratch = scratchDir_;
    if (scratch.empty()) {
        if (const char *env = std::getenv("_CONDOR_SCRATCH_DIR"))
            scratch = env;
        else
            scratch = ".";
    }

    char buf[4096];
    std::string absoluteScratch =
        realpath(scratch.c_str(), buf) ? std::string(buf) : scratch;

    auto perfProduct = std::make_unique<image::InferencePerfProduct>();

    for (auto const &cfg : models_) {
        std::string assets = cfg.assets.empty() ? assetsBaseDir_ : cfg.assets;
        std::string wrapper =
            cfg.wrapper.empty() ? defaultWrapper_ : cfg.wrapper;
        wrapper = resolve_under(assets, wrapper);
        std::string weights = resolve_under(assets, cfg.weights);

        mf::LogInfo("InferenceProduction")
            << "Running model: " << (cfg.label.empty() ? cfg.arch : cfg.label);

        auto result = image::InferenceProduction::runInferenceDetailed(
            detector_images, absoluteScratch, assets, cfg.arch, weights,
            wrapper, assets, cfg.want_cls, cfg.want_seg);

        image::ModelPerf perf;
        perf.model = cfg.label.empty() ? cfg.arch : cfg.label;
        perf.t_write_req_ms = static_cast<float>(result.perf.t_write_req_ms);
        perf.t_exec_total_ms = static_cast<float>(result.perf.t_exec_total_ms);
        perf.t_read_resp_ms = static_cast<float>(result.perf.t_read_resp_ms);
        perf.t_child_total_ms =
            static_cast<float>(result.perf.t_child_total_ms);
        perf.t_child_setup_ms =
            static_cast<float>(result.perf.t_child_setup_ms);
        perf.t_child_infer_ms =
            static_cast<float>(result.perf.t_child_infer_ms);
        perf.t_child_post_ms = static_cast<float>(result.perf.t_child_post_ms);
        perf.child_max_rss_mb =
            static_cast<float>(result.perf.child_max_rss_mb);
        perf.child_cuda_mem_mb =
            static_cast<float>(result.perf.child_cuda_mem_mb);
        perfProduct->per_model.push_back(std::move(perf));

    }

    e.put(std::move(perfProduct));
}

DEFINE_ART_MODULE(image::InferenceProducerModule)

} // namespace image
