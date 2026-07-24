#ifndef IMAGEPIPELINE_INFERENCEPROCESS_H
#define IMAGEPIPELINE_INFERENCEPROCESS_H

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>

namespace image {

/// Unique request/result/log paths retained on failure and removed on success.
class InferenceArtefacts {
  public:
    explicit InferenceArtefacts(const std::string &scratch_dir);
    ~InferenceArtefacts() noexcept;

    InferenceArtefacts(const InferenceArtefacts &) = delete;
    InferenceArtefacts &operator=(const InferenceArtefacts &) = delete;
    InferenceArtefacts(InferenceArtefacts &&) = delete;
    InferenceArtefacts &operator=(InferenceArtefacts &&) = delete;

    const std::string &requestPath() const noexcept { return request_path_; }
    const std::string &resultPath() const noexcept { return result_path_; }
    const std::string &metricsPath() const noexcept { return metrics_path_; }
    const std::string &stdoutPath() const noexcept { return stdout_path_; }
    const std::string &stderrPath() const noexcept { return stderr_path_; }

    std::string describe() const;
    void cleanupOnDestruction() noexcept { cleanup_on_destruction_ = true; }

  private:
    void cleanup() noexcept;

    std::string request_path_;
    std::string result_path_;
    std::string metrics_path_;
    std::string stdout_path_;
    std::string stderr_path_;
    bool cleanup_on_destruction_{false};
};

struct InferenceProcessConfig {
    std::string context;
    std::string arch;
    std::string weights_file;
    std::string inference_wrapper;
    std::string runtime_base_dir;
    std::string container_image;
    std::chrono::milliseconds timeout{0};
    std::size_t diagnostic_bytes{64U * 1024U};
    uint32_t width{0U};
    uint32_t height{0U};
};

struct InferenceExecution {
    double elapsed_ms{0.0};
};

/// Builds argv and runs the inference child directly with posix_spawnp.
class InferenceProcess {
  public:
    static std::string commandForLog(
        const InferenceProcessConfig &config,
        const InferenceArtefacts &artefacts);

    static InferenceExecution execute(const InferenceProcessConfig &config,
                                      const InferenceArtefacts &artefacts);

    static std::string diagnostics(const InferenceArtefacts &artefacts,
                                   std::size_t max_bytes);
};

} // namespace image

#endif // IMAGEPIPELINE_INFERENCEPROCESS_H
