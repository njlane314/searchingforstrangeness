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
#include "Products/InferencePred.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <algorithm>
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

std::string resolve_under(const std::string &base, const std::string &value) {
    if (value.empty() || base.empty())
        return value;

    if (value.front() == '/' ||
        (value.size() > 1 && value[1] == ':') ||
        value.find("://") != std::string::npos)
        return value;

    return join_path(base, value);
}

std::string resolve_scratch_dir(const std::string &configured) {
    std::string scratch = configured;
    if (scratch.empty()) {
        if (const char *env = std::getenv("_CONDOR_SCRATCH_DIR"))
            scratch = env;
        else
            scratch = ".";
    }

    char buf[4096];
    if (realpath(scratch.c_str(), buf))
        return std::string(buf);

    return scratch;
}

const image::ImageProduct *find_view(std::vector<image::ImageProduct> const &planes,
                                     geo::View_t view) {
    for (auto const &plane : planes) {
        if (plane.view == static_cast<int>(view))
            return &plane;
    }
    return nullptr;
}

} // namespace

class InferenceProducer : public art::EDProducer {
  public:
    struct ModelConfig {
        std::string label;
        std::string arch;
        std::string weights;
        std::string wrapper;
        std::string assets;
    };

    explicit InferenceProducer(fhicl::ParameterSet const &p);

    void produce(art::Event &e) override;

  private:
    static ModelConfig makeConfig(fhicl::ParameterSet const &p);

    art::InputTag planes_tag_;
    std::string scratch_dir_;
    std::string assets_base_dir_;
    std::string default_wrapper_;
    std::vector<ModelConfig> models_;
};

InferenceProducer::ModelConfig
InferenceProducer::makeConfig(fhicl::ParameterSet const &p) {
    ModelConfig cfg;
    cfg.label = p.get<std::string>("Label", "");
    cfg.arch = p.get<std::string>("Arch");
    cfg.weights = p.get<std::string>("Weights", "");
    cfg.wrapper = p.get<std::string>("Wrapper", "");
    cfg.assets = p.get<std::string>("Assets", "");
    return cfg;
}

InferenceProducer::InferenceProducer(fhicl::ParameterSet const &p)
    : planes_tag_{p.get<std::string>("PlanesTag")},
      scratch_dir_{p.get<std::string>("ScratchDir", "")},
      assets_base_dir_{p.get<std::string>("AssetsBaseDir", "")},
      default_wrapper_{p.get<std::string>("DefaultWrapper",
                                          "scripts/inference_wrapper.sh")} {
    produces<image::InferencePerfProduct>();
    produces<image::InferencePredProduct>();

    auto modelSets = p.get<std::vector<fhicl::ParameterSet>>("Models", {});
    models_.reserve(modelSets.size());
    for (auto const &modelPs : modelSets) {
        models_.push_back(makeConfig(modelPs));
    }
}

void InferenceProducer::produce(art::Event &e) {
    auto perfProduct = std::make_unique<image::InferencePerfProduct>();
    auto predProduct = std::make_unique<image::InferencePredProduct>();

    auto planes = e.getValidHandle<std::vector<image::ImageProduct>>(planes_tag_);
    if (planes->size() < 3) {
        // This happens, for example, when ImageProducer produced an empty
        // NuSlice (no neutrino slice). In that case we simply skip inference
        // and leave the perf/pred products empty.
        mf::LogDebug("InferenceProduction")
            << "Need at least three planes (U,V,W), got " << planes->size()
            << "; skipping inference for event " << e.id();

        e.put(std::move(perfProduct));
        e.put(std::move(predProduct));
        return;
    }

    auto const *U = find_view(*planes, geo::kU);
    auto const *V = find_view(*planes, geo::kV);
    auto const *W = find_view(*planes, geo::kW);
    if (!U || !V || !W) {
        throw cet::exception("InferenceProduction")
            << "Unable to identify U/V/W planes";
    }

    auto check_dims = [](image::ImageProduct const *p, char label) {
        std::size_t expected =
            static_cast<std::size_t>(p->width) *
            static_cast<std::size_t>(p->height);
        if (p->adc.size() != expected) {
            throw cet::exception("InferenceProduction")
                << "Plane " << label << " has incompatible dimensions: "
                << "width=" << p->width << " height=" << p->height
                << " adc.size()=" << p->adc.size()
                << " (expected " << expected << ")";
        }
    };

    check_dims(U, 'U');
    check_dims(V, 'V');
    check_dims(W, 'W');

    std::vector<image::ImageProduct> detector_images;
    detector_images.reserve(3);
    detector_images.emplace_back(*U);
    detector_images.emplace_back(*V);
    detector_images.emplace_back(*W);

    // Short-circuit if all detector images are effectively empty.
    // This avoids spinning up the external inference stack when
    // there is nothing in the images (e.g. no neutrino slice).
    auto is_empty = [](image::ImageProduct const &p) {
        if (p.adc.empty())
            return true;
        return std::all_of(p.adc.begin(), p.adc.end(),
                           [](float v) { return v == 0.0f; });
    };

    if (is_empty(detector_images[0]) &&
        is_empty(detector_images[1]) &&
        is_empty(detector_images[2])) {
        mf::LogDebug("InferenceProduction")
            << "Detector images are empty; skipping inference for event "
            << e.id();

        e.put(std::move(perfProduct));
        e.put(std::move(predProduct));
        return;
    }

    std::string absoluteScratch = resolve_scratch_dir(scratch_dir_);

    for (auto const &cfg : models_) {
        std::string assets = cfg.assets.empty() ? assets_base_dir_ : cfg.assets;
        std::string wrapper =
            cfg.wrapper.empty() ? default_wrapper_ : cfg.wrapper;
        wrapper = resolve_under(assets, wrapper);
        std::string weights = resolve_under(assets, cfg.weights);

        auto label_or_arch = cfg.label.empty() ? cfg.arch : cfg.label;
        std::string log_model = cfg.arch;
        if (!cfg.label.empty())
            log_model = cfg.label + " (arch: " + cfg.arch + ")";

        mf::LogInfo("InferenceProduction")
            << "Running model: " << log_model;

        auto result = image::InferenceProduction::runInference(
            detector_images, absoluteScratch, cfg.arch, weights, wrapper,
            assets);

        image::ModelPerf perf;
        perf.model = label_or_arch;
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
        perfProduct->per_model.push_back(std::move(perf));

        image::ModelPred pred;
        pred.model = label_or_arch;
        pred.scores = std::move(result.cls);
        predProduct->per_model.push_back(std::move(pred));

    }

    e.put(std::move(perfProduct));
    e.put(std::move(predProduct));
}

DEFINE_ART_MODULE(InferenceProducer)
