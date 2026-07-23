#include "ImagePipeline/InferenceProduction.h"

#include "ImagePipeline/InferenceProcess.h"
#include "ImagePipeline/InferenceProtocol.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <chrono>
#include <memory>
#include <stdexcept>

namespace image {
namespace {

[[noreturn]] void stageError(const InferenceProduction::Config &config,
                             const InferenceArtifacts &artifacts,
                             const std::string &stage,
                             const std::exception &error) {
    throw std::runtime_error(
        "Inference failed for " + config.context + " while " + stage + ": " +
        error.what() + "\n" +
        InferenceProcess::diagnostics(artifacts, config.diagnostic_bytes) +
        "\nArtifacts retained for debugging:\n" + artifacts.describe());
}

} // namespace

InferenceProduction::Result InferenceProduction::runInference(
    const ValidatedPlaneTriplet &detector_images, const Config &config) {
    std::unique_ptr<InferenceArtifacts> artifact_owner;
    try {
        artifact_owner =
            std::make_unique<InferenceArtifacts>(config.scratch_dir);
    } catch (const std::exception &error) {
        throw std::runtime_error(
            "Inference failed for " + config.context +
            " while creating temporary artifacts: " + error.what());
    }
    auto &artifacts = *artifact_owner;

    const auto write_start = std::chrono::steady_clock::now();
    try {
        InferenceProtocol::writeRequest(artifacts.requestPath(),
                                        detector_images);
    } catch (const std::exception &error) {
        stageError(config, artifacts, "writing the IASP request", error);
    }
    const auto write_end = std::chrono::steady_clock::now();

    InferenceProcessConfig process_config;
    process_config.context = config.context;
    process_config.arch = config.arch;
    process_config.weights_file = config.weights_file;
    process_config.inference_wrapper = config.inference_wrapper;
    process_config.runtime_base_dir = config.runtime_base_dir;
    process_config.container_image = config.container_image;
    process_config.timeout = config.timeout;
    process_config.diagnostic_bytes = config.diagnostic_bytes;
    process_config.width = detector_images.width();
    process_config.height = detector_images.height();

    try {
        mf::LogInfo("InferenceProduction")
            << "Executing inference for " << config.context << ": "
            << InferenceProcess::commandForLog(process_config, artifacts);
    } catch (const std::exception &error) {
        stageError(config, artifacts, "building the child command", error);
    }

    // execute() reports bounded stdout/stderr and retained paths on every
    // process-level failure.
    const auto execution =
        InferenceProcess::execute(process_config, artifacts);

    const auto read_start = std::chrono::steady_clock::now();
    Result result;
    InferenceChildMetrics child_metrics;
    try {
        result.cls = InferenceProtocol::readResult(
            artifacts.resultPath(), config.max_logits);
        child_metrics =
            InferenceProtocol::readMetrics(artifacts.metricsPath());
    } catch (const std::exception &error) {
        stageError(config, artifacts, "reading child output", error);
    }
    const auto read_end = std::chrono::steady_clock::now();

    result.perf.t_write_req_ms =
        std::chrono::duration<double, std::milli>(write_end - write_start)
            .count();
    result.perf.t_exec_total_ms = execution.elapsed_ms;
    result.perf.t_read_resp_ms =
        std::chrono::duration<double, std::milli>(read_end - read_start)
            .count();
    result.perf.t_child_total_ms = child_metrics.total_ms;
    result.perf.t_child_setup_ms = child_metrics.setup_ms;
    result.perf.t_child_infer_ms = child_metrics.infer_ms;
    result.perf.t_child_post_ms = child_metrics.post_ms;
    result.perf.child_max_rss_mb = child_metrics.max_rss_mb;

    mf::LogInfo("InferenceProduction")
        << "Inference complete for " << config.context
        << ": t_write_req_ms=" << result.perf.t_write_req_ms
        << " t_exec_total_ms=" << result.perf.t_exec_total_ms
        << " t_read_resp_ms=" << result.perf.t_read_resp_ms
        << " child(total/setup/infer/post)_ms="
        << result.perf.t_child_total_ms << "/"
        << result.perf.t_child_setup_ms << "/"
        << result.perf.t_child_infer_ms << "/"
        << result.perf.t_child_post_ms
        << " child_max_rss_mb=" << result.perf.child_max_rss_mb;

    if (!result.cls.empty()) {
        mf::LogInfo("InferenceProduction")
            << "First class score for " << config.context << ": "
            << result.cls.front();
    } else {
        mf::LogWarning("InferenceProduction")
            << "No class scores returned for " << config.context;
    }

    artifacts.cleanupOnDestruction();
    return result;
}

} // namespace image
