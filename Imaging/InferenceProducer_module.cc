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
#include "Products/InferenceResult.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

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
    };

    explicit InferenceProducerModule(fhicl::ParameterSet const &p);

    void produce(art::Event &e) override;

  private:
    static ModelConfig make_config(fhicl::ParameterSet const &p);

    art::InputTag planes_tag_;
    std::string scratch_dir_;
    std::string assets_base_dir_;
    std::string default_wrapper_;
    std::vector<ModelConfig> models_;
};

InferenceProducerModule::ModelConfig
InferenceProducerModule::make_config(fhicl::ParameterSet const &p) {
    ModelConfig model_config;
    model_config.label = p.get<std::string>("Label", "");
    model_config.arch = p.get<std::string>("Arch");
    model_config.weights = p.get<std::string>("Weights", "");
    model_config.wrapper = p.get<std::string>("Wrapper", "");
    model_config.assets = p.get<std::string>("Assets", "");
    return model_config;
}

InferenceProducerModule::InferenceProducerModule(
    fhicl::ParameterSet const &p)
    : planes_tag_{p.get<std::string>("PlanesTag")},
      scratch_dir_{p.get<std::string>("ScratchDir", "")},
      assets_base_dir_{p.get<std::string>("AssetsBaseDir", "")},
      default_wrapper_{resolve_under(
          assets_base_dir_, p.get<std::string>("DefaultWrapper",
                                              "scripts/inference_wrapper.sh"))} {
    produces<image::InferencePerfProduct>();
    produces<image::InferenceResultProduct>();

    auto model_parameter_sets =
        p.get<std::vector<fhicl::ParameterSet>>("Models", {});
    models_.reserve(model_parameter_sets.size());
    for (auto const &model_parameter_set : model_parameter_sets) {
        models_.push_back(make_config(model_parameter_set));
    }
}

void InferenceProducerModule::produce(art::Event &e) {
    auto plane_handle =
        e.getValidHandle<std::vector<image::PlaneImage>>(planes_tag_);
    if (plane_handle->size() < 3) {
        throw cet::exception("InferenceProduction")
            << "Need at least three planes (U,V,W), got "
            << plane_handle->size();
    }

    auto const *first_plane = &plane_handle->front();
    auto const *plane_u = find_view(*plane_handle, geo::kU, first_plane);
    auto const *plane_v = find_view(*plane_handle, geo::kV, first_plane);
    auto const *plane_w = find_view(*plane_handle, geo::kW, first_plane);
    if (plane_u == nullptr || plane_v == nullptr || plane_w == nullptr) {
        throw cet::exception("InferenceProduction")
            << "Unable to identify U/V/W planes";
    }

    image::PlaneImage plane_image_u = *plane_u;
    image::PlaneImage plane_image_v = *plane_v;
    image::PlaneImage plane_image_w = *plane_w;

    std::vector<image::PlaneImage> detector_images;
    detector_images.reserve(3);
    detector_images.emplace_back(std::move(plane_image_u));
    detector_images.emplace_back(std::move(plane_image_v));
    detector_images.emplace_back(std::move(plane_image_w));

    std::string scratch_directory = scratch_dir_;
    if (scratch_directory.empty()) {
        if (const char *scratch_env = std::getenv("_CONDOR_SCRATCH_DIR"))
            scratch_directory = scratch_env;
        else
            scratch_directory = ".";
    }

    char scratch_buffer[4096];
    std::string absolute_scratch =
        realpath(scratch_directory.c_str(), scratch_buffer)
            ? std::string(scratch_buffer)
            : scratch_directory;

    auto perf_product = std::make_unique<image::InferencePerfProduct>();
    auto result_product = std::make_unique<image::InferenceResultProduct>();

    for (auto const &model_config : models_) {
        std::string assets_dir =
            model_config.assets.empty() ? assets_base_dir_ : model_config.assets;
        std::string wrapper_path =
            model_config.wrapper.empty() ? default_wrapper_ : model_config.wrapper;
        wrapper_path = resolve_under(assets_dir, wrapper_path);
        std::string weights_path = resolve_under(assets_dir, model_config.weights);

        mf::LogInfo("InferenceProduction")
            << "Running model: "
            << (model_config.label.empty() ? model_config.arch
                                           : model_config.label);

        auto inference_result = image::InferenceProduction::runInferenceDetailed(
            detector_images, absolute_scratch, assets_dir, model_config.arch,
            weights_path, wrapper_path, assets_dir);

        image::ModelPerf model_perf;
        model_perf.model =
            model_config.label.empty() ? model_config.arch : model_config.label;
        model_perf.t_write_req_ms =
            static_cast<float>(inference_result.perf.t_write_req_ms);
        model_perf.t_exec_total_ms =
            static_cast<float>(inference_result.perf.t_exec_total_ms);
        model_perf.t_read_resp_ms =
            static_cast<float>(inference_result.perf.t_read_resp_ms);
        model_perf.t_child_total_ms =
            static_cast<float>(inference_result.perf.t_child_total_ms);
        model_perf.t_child_setup_ms =
            static_cast<float>(inference_result.perf.t_child_setup_ms);
        model_perf.t_child_infer_ms =
            static_cast<float>(inference_result.perf.t_child_infer_ms);
        model_perf.t_child_post_ms =
            static_cast<float>(inference_result.perf.t_child_post_ms);
        model_perf.child_max_rss_mb =
            static_cast<float>(inference_result.perf.child_max_rss_mb);
        perf_product->per_model.push_back(model_perf);

        image::ModelResult model_result;
        model_result.model =
            model_config.label.empty() ? model_config.arch : model_config.label;
        model_result.cls_scores = std::move(inference_result.cls);
        model_result.segmentation_width = inference_result.segW;
        model_result.segmentation_height = inference_result.segH;
        model_result.segmentation_u = std::move(inference_result.seg_u);
        model_result.segmentation_v = std::move(inference_result.seg_v);
        model_result.segmentation_w = std::move(inference_result.seg_w);
        model_result.confidence_u = std::move(inference_result.conf_u);
        model_result.confidence_v = std::move(inference_result.conf_v);
        model_result.confidence_w = std::move(inference_result.conf_w);
        model_result.features = std::move(inference_result.features);
        model_result.feature_seed = inference_result.feature_seed;
        result_product->per_model.push_back(std::move(model_result));

    }

    e.put(std::move(perf_product));
    e.put(std::move(result_product));
}

DEFINE_ART_MODULE(image::InferenceProducerModule)

} // namespace image
