#include "ImagePipeline/InferenceConfiguration.h"

#include "cetlib_except/exception.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <utility>

namespace image {
namespace {

constexpr const char *kDefaultContainer =
    "/cvmfs/uboone.opensciencegrid.org/containers/"
    "lantern_v2_me_06_03_prod";
constexpr uint32_t kDefaultDiagnosticBytes = 64U * 1024U;
constexpr uint32_t kMaximumDiagnosticBytes = 4U * 1024U * 1024U;
constexpr uint32_t kDefaultMaxLogits = 1024U * 1024U;
constexpr uint32_t kMaximumMaxLogits = 16U * 1024U * 1024U;

std::string joinPath(std::string base, const std::string &leaf) {
    if (leaf.empty()) {
        return base;
    }
    if (base.empty()) {
        return leaf;
    }
    if (base.back() != '/') {
        base.push_back('/');
    }
    return base + leaf;
}

std::string resolveUnder(const std::string &base,
                         const std::string &value) {
    if (value.empty() || base.empty()) {
        return value;
    }
    if (value.front() == '/' ||
        (value.size() > 1U && value[1] == ':') ||
        value.find("://") != std::string::npos) {
        return value;
    }
    return joinPath(base, value);
}

std::string resolveScratchDirectory(const std::string &configured) {
    std::string scratch = configured;
    if (scratch.empty()) {
        if (const char *condor_scratch = std::getenv("_CONDOR_SCRATCH_DIR")) {
            scratch = condor_scratch;
        } else {
            scratch = ".";
        }
    }

    char resolved[4096];
    return ::realpath(scratch.c_str(), resolved) != nullptr
               ? std::string(resolved)
               : scratch;
}

template <typename T>
T getFirstPresent(const fhicl::ParameterSet &parameters,
                  std::initializer_list<const char *> keys, T fallback) {
    for (const auto *key : keys) {
        if (parameters.has_key(key)) {
            return parameters.get<T>(key);
        }
    }
    return fallback;
}

std::chrono::milliseconds timeoutFromSeconds(const double seconds) {
    if (!std::isfinite(seconds) || seconds < 0.0) {
        throw cet::exception("InferenceProduction")
            << "TimeoutSeconds must be finite and non-negative";
    }
    if (seconds == 0.0) {
        return std::chrono::milliseconds{0};
    }

    const double milliseconds = std::ceil(seconds * 1000.0);
    if (milliseconds >
        static_cast<double>(
            std::numeric_limits<std::chrono::milliseconds::rep>::max())) {
        throw cet::exception("InferenceProduction")
            << "TimeoutSeconds is too large";
    }
    return std::chrono::milliseconds{
        static_cast<std::chrono::milliseconds::rep>(milliseconds)};
}

std::string resolveWeightsPath(const std::string &runtime_dir,
                               const InferenceModelConfig &model) {
    if (!model.weights.empty()) {
        return resolveUnder(runtime_dir, model.weights);
    }
    if (!model.checkpoint.empty()) {
        return resolveUnder(
            runtime_dir,
            joinPath(model.weights_dir, model.checkpoint));
    }
    return {};
}

void validateModel(const InferenceModelConfig &model) {
    if (model.arch.empty()) {
        throw cet::exception("InferenceProduction")
            << "Each model requires Arch or ModelDefaults.Arch";
    }
    if (model.weights.empty() && model.checkpoint.empty()) {
        throw cet::exception("InferenceProduction")
            << "Each model requires Weights or Checkpoint";
    }
    if (model.wrapper.empty()) {
        throw cet::exception("InferenceProduction")
            << "Each model requires an inference Wrapper";
    }
    if (model.container.empty()) {
        throw cet::exception("InferenceProduction")
            << "Each model requires a Container image";
    }
    timeoutFromSeconds(model.timeout_seconds);
    if (model.diagnostic_bytes == 0U ||
        model.diagnostic_bytes > kMaximumDiagnosticBytes) {
        throw cet::exception("InferenceProduction")
            << "DiagnosticBytes must be in [1, "
            << kMaximumDiagnosticBytes << "]";
    }
    if (model.max_logits == 0U ||
        model.max_logits > kMaximumMaxLogits) {
        throw cet::exception("InferenceProduction")
            << "MaxLogits must be in [1, " << kMaximumMaxLogits << "]";
    }
}

} // namespace

std::string InferenceModelConfig::productName() const {
    return label.empty() ? arch : label;
}

InferenceModelConfig InferenceConfiguration::mergeModelConfig(
    const fhicl::ParameterSet &parameters,
    InferenceModelConfig defaults) {
    InferenceModelConfig model = std::move(defaults);
    const bool has_weights = parameters.has_key("Weights");
    const bool has_weights_dir = parameters.has_key("WeightsDir");
    const bool has_checkpoint = parameters.has_key("Checkpoint");

    model.label =
        getFirstPresent<std::string>(parameters, {"Label"}, model.label);
    model.arch =
        getFirstPresent<std::string>(parameters, {"Arch"}, model.arch);
    model.weights =
        getFirstPresent<std::string>(parameters, {"Weights"}, model.weights);
    model.weights_dir = getFirstPresent<std::string>(
        parameters, {"WeightsDir"}, model.weights_dir);
    model.checkpoint = getFirstPresent<std::string>(
        parameters, {"Checkpoint"}, model.checkpoint);
    model.wrapper =
        getFirstPresent<std::string>(parameters, {"Wrapper"}, model.wrapper);
    model.runtime =
        getFirstPresent<std::string>(parameters, {"Runtime"}, model.runtime);
    model.container = getFirstPresent<std::string>(
        parameters, {"Container", "ContainerImage"}, model.container);
    model.timeout_seconds = getFirstPresent<double>(
        parameters, {"TimeoutSeconds"}, model.timeout_seconds);
    model.diagnostic_bytes = getFirstPresent<uint32_t>(
        parameters, {"DiagnosticBytes"}, model.diagnostic_bytes);
    model.max_logits = getFirstPresent<uint32_t>(
        parameters, {"MaxLogits"}, model.max_logits);

    if (!has_weights && (has_weights_dir || has_checkpoint)) {
        model.weights.clear();
    }
    return model;
}

InferenceConfiguration::InferenceConfiguration(
    const fhicl::ParameterSet &parameters)
    : runtime_base_dir_{
          parameters.get<std::string>("RuntimeBaseDir", "")},
      default_wrapper_{parameters.get<std::string>(
          "DefaultInferenceWrapper", "scripts/inference_wrapper.sh")},
      scratch_dir_{resolveScratchDirectory(
          parameters.get<std::string>("ScratchDir", ""))} {
    InferenceModelConfig defaults;
    defaults.wrapper = default_wrapper_;
    defaults.runtime = runtime_base_dir_;
    defaults.container =
        parameters.get<std::string>("ContainerImage", kDefaultContainer);
    defaults.timeout_seconds =
        parameters.get<double>("TimeoutSeconds", 0.0);
    defaults.diagnostic_bytes = parameters.get<uint32_t>(
        "DiagnosticBytes", kDefaultDiagnosticBytes);
    defaults.max_logits =
        parameters.get<uint32_t>("MaxLogits", kDefaultMaxLogits);
    if (parameters.has_key("ModelDefaults")) {
        defaults = mergeModelConfig(
            parameters.get<fhicl::ParameterSet>("ModelDefaults"),
            std::move(defaults));
    }

    const auto model_sets =
        parameters.get<std::vector<fhicl::ParameterSet>>("Models", {});
    if (model_sets.empty()) {
        throw cet::exception("InferenceProduction")
            << "InferenceProducer requires at least one Models entry";
    }

    models_.reserve(model_sets.size());
    for (const auto &model_parameters : model_sets) {
        auto model = mergeModelConfig(model_parameters, defaults);
        validateModel(model);
        models_.push_back(std::move(model));
    }
}

InferenceProduction::Config InferenceConfiguration::makeRunConfig(
    const InferenceModelConfig &model,
    const std::string &event_context) const {
    const std::string runtime =
        model.runtime.empty() ? runtime_base_dir_ : model.runtime;
    const std::string wrapper = resolveUnder(
        runtime, model.wrapper.empty() ? default_wrapper_ : model.wrapper);

    std::ostringstream context;
    context << event_context << ", model " << model.productName()
            << " (arch: " << model.arch << ")";

    InferenceProduction::Config config;
    config.scratch_dir = scratch_dir_;
    config.context = context.str();
    config.arch = model.arch;
    config.weights_file = resolveWeightsPath(runtime, model);
    config.inference_wrapper = wrapper;
    config.runtime_base_dir = runtime;
    config.container_image = model.container;
    config.timeout = timeoutFromSeconds(model.timeout_seconds);
    config.diagnostic_bytes = model.diagnostic_bytes;
    config.max_logits = model.max_logits;
    return config;
}

} // namespace image
