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
                             const InferenceArtefacts &artefacts,
                             const std::string &stage,
                             const std::exception &error) {
    throw std::runtime_error(
        "Inference failed for " + config.context + " while " + stage + ": " +
        error.what() + "\n" +
        InferenceProcess::diagnostics(artefacts, config.diagnostic_bytes) +
        "\nArtefacts retained for debugging:\n" + artefacts.describe());
}

} // namespace

InferenceProduction::Result InferenceProduction::runInference(
    const ValidatedPlaneTriplet &detector_images, const Config &config) {
    std::unique_ptr<InferenceArtefacts> artefact_owner;
    try {
        artefact_owner =
            std::make_unique<InferenceArtefacts>(config.scratch_dir);
    } catch (const std::exception &error) {
        throw std::runtime_error(
            "Inference failed for " + config.context +
            " while creating temporary artefacts: " + error.what());
    }
    auto &artefacts = *artefact_owner;

    const auto write_start = std::chrono::steady_clock::now();
    try {
        InferenceProtocol::writeRequest(artefacts.requestPath(),
                                        detector_images);
    } catch (const std::exception &error) {
        stageError(config, artefacts, "writing the IASP request", error);
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
            << InferenceProcess::commandForLog(process_config, artefacts);
    } catch (const std::exception &error) {
        stageError(config, artefacts, "building the child command", error);
    }

    // execute() reports bounded stdout/stderr and retained paths on every
    // process-level failure.
    const auto execution =
        InferenceProcess::execute(process_config, artefacts);

    const auto read_start = std::chrono::steady_clock::now();
    Result result;
    InferenceChildMetrics child_metrics;
    try {
        result.cls = InferenceProtocol::readResult(
            artefacts.resultPath(), config.max_logits);
        child_metrics =
            InferenceProtocol::readMetrics(artefacts.metricsPath());
    } catch (const std::exception &error) {
        stageError(config, artefacts, "reading child output", error);
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

    artefacts.cleanupOnDestruction();
    return result;
}

} // namespace image
