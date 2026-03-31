#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "cetlib_except/exception.h"

#include "ImagePipeline/InferenceProduction.h"
#include "Products/SparsePlaneImage.h"
#include "Products/InferenceMetrics.h"
#include "Products/InferencePredictions.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
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

const image::SparsePlaneImage *find_view(
    std::vector<image::SparsePlaneImage> const &planes, geo::View_t view) {
    for (auto const &plane : planes) {
        if (plane.view == static_cast<int>(view))
            return &plane;
    }
    return nullptr;
}

}

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
    produces<image::InferenceMetrics>();
    produces<image::InferencePredictions>();

    auto modelSets = p.get<std::vector<fhicl::ParameterSet>>("Models", {});
    models_.reserve(modelSets.size());
    for (auto const &modelPs : modelSets) {
        models_.push_back(makeConfig(modelPs));
    }
}

void InferenceProducer::produce(art::Event &e) {
    auto perfProduct = std::make_unique<image::InferenceMetrics>();
    auto predProduct = std::make_unique<image::InferencePredictions>();

    auto planes = e.getValidHandle<std::vector<image::SparsePlaneImage>>(planes_tag_);
    if (planes->size() < 3) {
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

    auto check_plane = [](image::SparsePlaneImage const *p, char label) {
        if (p->feature_dim == 0) {
            throw cet::exception("InferenceProduction")
                << "Plane " << label << " has invalid feature_dim=0";
        }
        if ((p->coords.size() % 2) != 0) {
            throw cet::exception("InferenceProduction")
                << "Plane " << label << " has odd coords.size()=" << p->coords.size();
        }

        std::size_t const nnz = p->coords.size() / 2;
        if (p->features.size() != nnz * static_cast<std::size_t>(p->feature_dim)) {
            throw cet::exception("InferenceProduction")
                << "Plane " << label << " has mismatched sparse payload: "
                << "coords.size()=" << p->coords.size()
                << " feature_dim=" << p->feature_dim
                << " features.size()=" << p->features.size();
        }
        if (!p->semantic.empty() && p->semantic.size() != nnz) {
            throw cet::exception("InferenceProduction")
                << "Plane " << label << " has semantic.size()=" << p->semantic.size()
                << " but nnz=" << nnz;
        }
        for (std::size_t i = 0; i < nnz; ++i) {
            auto const row = p->coords[2 * i];
            auto const col = p->coords[2 * i + 1];
            if (row < 0 || col < 0 ||
                static_cast<std::size_t>(row) >= p->height ||
                static_cast<std::size_t>(col) >= p->width) {
                throw cet::exception("InferenceProduction")
                    << "Plane " << label << " has out-of-range sparse coordinate ("
                    << row << ", " << col << ") for image size "
                    << p->width << "x" << p->height;
            }
        }
    };

    check_plane(U, 'U');
    check_plane(V, 'V');
    check_plane(W, 'W');

    if (U->width != V->width || U->width != W->width ||
        U->height != V->height || U->height != W->height) {
        throw cet::exception("InferenceProduction")
            << "U/V/W planes must have matching dimensions, got "
            << "U=" << U->width << "x" << U->height << " "
            << "V=" << V->width << "x" << V->height << " "
            << "W=" << W->width << "x" << W->height;
    }

    std::vector<image::SparsePlaneImage> detector_images;
    detector_images.reserve(3);
    detector_images.emplace_back(*U);
    detector_images.emplace_back(*V);
    detector_images.emplace_back(*W);

    auto is_empty = [](image::SparsePlaneImage const &p) {
        return p.coords.empty();
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

        image::ModelMetrics perf;
        perf.model = label_or_arch;
        perf.t_exec_total_ms = static_cast<float>(result.perf.t_exec_total_ms);
        perf.t_child_infer_ms =
            static_cast<float>(result.perf.t_child_infer_ms);
        perf.child_max_rss_mb =
            static_cast<float>(result.perf.child_max_rss_mb);
        perfProduct->per_model.push_back(std::move(perf));

        image::ModelPrediction pred;
        pred.model = label_or_arch;
        pred.scores = std::move(result.cls);
        predProduct->per_model.push_back(std::move(pred));

    }

    e.put(std::move(perfProduct));
    e.put(std::move(predProduct));
}

DEFINE_ART_MODULE(InferenceProducer)
