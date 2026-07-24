#include "ImagePipeline/InferenceProcess.h"

#include <algorithm>
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <spawn.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

extern char **environ;

namespace image {
namespace {

constexpr const char *kApptainerExecutable = "apptainer";
constexpr std::chrono::milliseconds kWaitPollInterval{50};
constexpr std::chrono::seconds kTerminateGracePeriod{2};

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

bool isAbsoluteOrUri(const std::string &path) {
    return !path.empty() &&
           (path.front() == '/' ||
            (path.size() > 1U && path[1] == ':') ||
            path.find("://") != std::string::npos);
}

std::string makeAbsolutePath(const std::string &path) {
    if (path.empty() || isAbsoluteOrUri(path)) {
        return path;
    }

    char resolved[4096];
    if (::realpath(path.c_str(), resolved) != nullptr) {
        return std::string(resolved);
    }

    char cwd[4096];
    if (::getcwd(cwd, sizeof(cwd)) != nullptr) {
        return joinPath(std::string(cwd), path);
    }
    return path;
}

bool isLocalFilesystemPath(const std::string &path) {
    return !path.empty() && path.front() == '/' &&
           path.find("://") == std::string::npos;
}

std::string parentDirectory(const std::string &path) {
    if (!isLocalFilesystemPath(path)) {
        return {};
    }
    const auto separator = path.rfind('/');
    if (separator == std::string::npos) {
        return {};
    }
    return separator == 0U ? "/" : path.substr(0U, separator);
}

void appendUnique(std::vector<std::string> &items,
                  std::set<std::string> &seen, const std::string &item) {
    if (!item.empty() && seen.insert(item).second) {
        if (item.find(',') != std::string::npos) {
            throw std::invalid_argument(
                "Apptainer bind paths cannot contain commas: " + item);
        }
        items.push_back(item);
    }
}

bool cvmfsAlreadyAvailable() {
    if (const char *bind_path = std::getenv("APPTAINER_BINDPATH")) {
        return std::string(bind_path).find("/cvmfs") != std::string::npos;
    }

    std::ifstream mounts("/proc/mounts");
    std::string line;
    while (std::getline(mounts, line)) {
        if (line.find(" /cvmfs ") != std::string::npos) {
            return true;
        }
    }
    return false;
}

void validateConfig(const InferenceProcessConfig &config) {
    if (config.context.empty()) {
        throw std::invalid_argument("Inference process context is empty");
    }
    if (config.arch.empty()) {
        throw std::invalid_argument("Inference architecture is empty");
    }
    if (config.weights_file.empty()) {
        throw std::invalid_argument("Inference weights path is empty");
    }
    if (config.inference_wrapper.empty()) {
        throw std::invalid_argument("Inference wrapper path is empty");
    }
    if (config.container_image.empty()) {
        throw std::invalid_argument("Inference container image is empty");
    }
    if (config.width == 0U || config.height == 0U) {
        throw std::invalid_argument(
            "Inference image dimensions must be non-zero");
    }
    if (config.timeout.count() < 0) {
        throw std::invalid_argument(
            "Inference process timeout cannot be negative");
    }
    if (config.diagnostic_bytes == 0U) {
        throw std::invalid_argument(
            "Inference diagnostic byte limit must be non-zero");
    }
}

std::vector<std::string>
buildArguments(const InferenceProcessConfig &config,
               const InferenceArtefacts &artefacts) {
    validateConfig(config);

    std::string runtime_dir = config.runtime_base_dir;
    if (runtime_dir.empty()) {
        const auto scripts = config.inference_wrapper.rfind("/scripts/");
        runtime_dir =
            scripts == std::string::npos
                ? std::string{}
                : config.inference_wrapper.substr(0U, scripts);
    }

    const std::string absolute_runtime = makeAbsolutePath(runtime_dir);
    const std::string absolute_wrapper =
        makeAbsolutePath(config.inference_wrapper);
    const std::string absolute_weights =
        makeAbsolutePath(config.weights_file);

    std::vector<std::string> binds;
    std::set<std::string> seen;
    appendUnique(binds, seen, parentDirectory(artefacts.requestPath()));
    appendUnique(binds, seen, absolute_runtime);
    appendUnique(binds, seen, parentDirectory(absolute_wrapper));
    appendUnique(binds, seen, parentDirectory(absolute_weights));
    if (!cvmfsAlreadyAvailable()) {
        binds.insert(binds.begin(), "/cvmfs");
    }

    std::ostringstream bind_csv;
    for (std::size_t index = 0; index < binds.size(); ++index) {
        if (index != 0U) {
            bind_csv << ',';
        }
        bind_csv << binds[index];
    }

    return {
        kApptainerExecutable,
        "exec",
        "--cleanenv",
        "--bind",
        bind_csv.str(),
        config.container_image,
        "/bin/bash",
        absolute_wrapper,
        "--in",
        artefacts.requestPath(),
        "--out",
        artefacts.resultPath(),
        "--metrics",
        artefacts.metricsPath(),
        "--W",
        std::to_string(config.width),
        "--H",
        std::to_string(config.height),
        "--arch",
        config.arch,
        "--weights",
        absolute_weights,
    };
}

std::string quoteForLog(const std::string &argument) {
    if (argument.find_first_of(" \t\n'\"") == std::string::npos) {
        return argument;
    }

    std::string quoted{"'"};
    for (const char character : argument) {
        if (character == '\'') {
            quoted += "'\\''";
        } else {
            quoted.push_back(character);
        }
    }
    quoted.push_back('\'');
    return quoted;
}

std::string joinArgumentsForLog(const std::vector<std::string> &arguments) {
    std::ostringstream output;
    for (std::size_t index = 0; index < arguments.size(); ++index) {
        if (index != 0U) {
            output << ' ';
        }
        output << quoteForLog(arguments[index]);
    }
    return output.str();
}

bool pathExists(const std::string &path) {
    struct stat info {};
    return ::stat(path.c_str(), &info) == 0;
}

std::string readTail(const std::string &path, const std::size_t max_bytes) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) {
        return "<not created>";
    }

    stream.seekg(0, std::ios::end);
    const std::streamoff size = stream.tellg();
    if (size < 0) {
        return "<could not determine capture size>";
    }

    const std::streamoff retained =
        std::min<std::streamoff>(size,
                                 static_cast<std::streamoff>(max_bytes));
    stream.seekg(size - retained, std::ios::beg);
    std::string content(static_cast<std::size_t>(retained), '\0');
    if (retained != 0) {
        stream.read(&content[0], retained);
        if (stream.gcount() != retained) {
            return "<could not read capture>";
        }
    }

    if (size > retained) {
        return "<truncated; showing last " + std::to_string(retained) +
               " bytes>\n" + content;
    }
    return content.empty() ? "<empty>" : content;
}

class SpawnFileActions {
  public:
    SpawnFileActions() {
        const int error = ::posix_spawn_file_actions_init(&actions_);
        if (error != 0) {
            throw std::runtime_error(
                "posix_spawn_file_actions_init failed: " +
                std::string(std::strerror(error)));
        }
        ready_ = true;
    }

    ~SpawnFileActions() noexcept {
        if (ready_) {
            ::posix_spawn_file_actions_destroy(&actions_);
        }
    }

    SpawnFileActions(const SpawnFileActions &) = delete;
    SpawnFileActions &operator=(const SpawnFileActions &) = delete;

    void redirect(const int descriptor, const std::string &path) {
        const int error = ::posix_spawn_file_actions_addopen(
            &actions_, descriptor, path.c_str(),
            O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (error != 0) {
            throw std::runtime_error("Could not redirect fd " +
                                     std::to_string(descriptor) + " to " +
                                     path + ": " + std::strerror(error));
        }
    }

    posix_spawn_file_actions_t *get() noexcept { return &actions_; }

  private:
    posix_spawn_file_actions_t actions_{};
    bool ready_{false};
};

class SpawnAttributes {
  public:
    SpawnAttributes() {
        int error = ::posix_spawnattr_init(&attributes_);
        if (error != 0) {
            throw std::runtime_error("posix_spawnattr_init failed: " +
                                     std::string(std::strerror(error)));
        }
        ready_ = true;

        error = ::posix_spawnattr_setpgroup(&attributes_, 0);
        if (error == 0) {
            short flags = 0;
            error = ::posix_spawnattr_getflags(&attributes_, &flags);
            if (error == 0) {
                error = ::posix_spawnattr_setflags(
                    &attributes_,
                    static_cast<short>(flags | POSIX_SPAWN_SETPGROUP));
            }
        }
        if (error != 0) {
            ::posix_spawnattr_destroy(&attributes_);
            ready_ = false;
            throw std::runtime_error(
                "Could not configure inference process group: " +
                std::string(std::strerror(error)));
        }
    }

    ~SpawnAttributes() noexcept {
        if (ready_) {
            ::posix_spawnattr_destroy(&attributes_);
        }
    }

    SpawnAttributes(const SpawnAttributes &) = delete;
    SpawnAttributes &operator=(const SpawnAttributes &) = delete;

    posix_spawnattr_t *get() noexcept { return &attributes_; }

  private:
    posix_spawnattr_t attributes_{};
    bool ready_{false};
};

void signalProcessGroup(const pid_t child, const int signal) noexcept {
    if (::kill(-child, signal) != 0 && errno != ESRCH) {
        // Best effort. The subsequent wait still determines child state.
    }
}

pid_t waitpidRetry(const pid_t child, int *status, const int options) {
    while (true) {
        const pid_t result = ::waitpid(child, status, options);
        if (result < 0 && errno == EINTR) {
            continue;
        }
        return result;
    }
}

class ChildGuard {
  public:
    explicit ChildGuard(const pid_t child) : child_{child} {}

    ~ChildGuard() noexcept {
        if (armed_) {
            signalProcessGroup(child_, SIGKILL);
            int status = 0;
            waitpidRetry(child_, &status, 0);
        }
    }

    void release() noexcept { armed_ = false; }

  private:
    pid_t child_{0};
    bool armed_{true};
};

struct WaitResult {
    int status{0};
    bool timed_out{false};
};

WaitResult waitForChild(const pid_t child,
                        const std::chrono::milliseconds timeout) {
    WaitResult result;
    if (timeout.count() == 0) {
        if (waitpidRetry(child, &result.status, 0) != child) {
            throw std::runtime_error("waitpid failed: " +
                                     std::string(std::strerror(errno)));
        }
        return result;
    }

    const auto deadline = std::chrono::steady_clock::now() + timeout;
    while (std::chrono::steady_clock::now() < deadline) {
        const pid_t waited = waitpidRetry(child, &result.status, WNOHANG);
        if (waited == child) {
            return result;
        }
        if (waited < 0) {
            throw std::runtime_error("waitpid failed: " +
                                     std::string(std::strerror(errno)));
        }
        std::this_thread::sleep_for(kWaitPollInterval);
    }

    result.timed_out = true;
    signalProcessGroup(child, SIGTERM);
    const auto terminate_deadline =
        std::chrono::steady_clock::now() + kTerminateGracePeriod;
    while (std::chrono::steady_clock::now() < terminate_deadline) {
        const pid_t waited = waitpidRetry(child, &result.status, WNOHANG);
        if (waited == child) {
            return result;
        }
        if (waited < 0) {
            throw std::runtime_error(
                "waitpid failed while terminating timed-out child: " +
                std::string(std::strerror(errno)));
        }
        std::this_thread::sleep_for(kWaitPollInterval);
    }

    signalProcessGroup(child, SIGKILL);
    if (waitpidRetry(child, &result.status, 0) != child) {
        throw std::runtime_error(
            "waitpid failed after killing timed-out child: " +
            std::string(std::strerror(errno)));
    }
    return result;
}

std::string describeStatus(const int status) {
    if (WIFEXITED(status)) {
        return "exit code " + std::to_string(WEXITSTATUS(status));
    }
    if (WIFSIGNALED(status)) {
        return "signal " + std::to_string(WTERMSIG(status));
    }
    return "wait status " + std::to_string(status);
}

[[noreturn]] void processError(const InferenceProcessConfig &config,
                               const InferenceArtefacts &artefacts,
                               const std::string &message) {
    throw std::runtime_error(
        "Inference process failed for " + config.context + ": " + message +
        "\n" +
        InferenceProcess::diagnostics(artefacts, config.diagnostic_bytes) +
        "\nArtefacts retained for debugging:\n" + artefacts.describe());
}

} // namespace

InferenceArtefacts::InferenceArtefacts(const std::string &scratch_dir) {
    const std::string directory = scratch_dir.empty() ? "." : scratch_dir;
    std::string pattern = joinPath(directory, "ia_XXXXXX");
    std::vector<char> writable(pattern.begin(), pattern.end());
    writable.push_back('\0');

    const int descriptor = ::mkstemp(writable.data());
    if (descriptor < 0) {
        throw std::runtime_error("Could not create inference artefacts in " +
                                 directory + ": " + std::strerror(errno));
    }
    if (::close(descriptor) != 0 && errno != EINTR) {
        const int saved_errno = errno;
        ::remove(writable.data());
        throw std::runtime_error(
            "Could not close inference request file " +
            std::string(writable.data()) + ": " +
            std::strerror(saved_errno));
    }

    request_path_ = writable.data();
    result_path_ = request_path_ + "_results.bin";
    metrics_path_ = result_path_ + ".meta";
    stdout_path_ = request_path_ + "_py_script.out";
    stderr_path_ = request_path_ + "_py_script.err";
}

InferenceArtefacts::~InferenceArtefacts() noexcept {
    if (cleanup_on_destruction_) {
        cleanup();
    }
}

void InferenceArtefacts::cleanup() noexcept {
    ::remove(request_path_.c_str());
    ::remove(result_path_.c_str());
    ::remove(metrics_path_.c_str());
    ::remove(stdout_path_.c_str());
    ::remove(stderr_path_.c_str());
}

std::string InferenceArtefacts::describe() const {
    std::ostringstream output;
    output << "  request: " << request_path_ << '\n'
           << "  result: " << result_path_ << '\n'
           << "  metrics: " << metrics_path_ << '\n'
           << "  stdout: " << stdout_path_ << '\n'
           << "  stderr: " << stderr_path_;
    return output.str();
}

std::string InferenceProcess::commandForLog(
    const InferenceProcessConfig &config,
    const InferenceArtefacts &artefacts) {
    return joinArgumentsForLog(buildArguments(config, artefacts));
}

InferenceExecution
InferenceProcess::execute(const InferenceProcessConfig &config,
                          const InferenceArtefacts &artefacts) {
    try {
        std::vector<std::string> arguments =
            buildArguments(config, artefacts);

        SpawnFileActions file_actions;
        SpawnAttributes attributes;
        file_actions.redirect(STDOUT_FILENO, artefacts.stdoutPath());
        file_actions.redirect(STDERR_FILENO, artefacts.stderrPath());

        std::vector<char *> argv;
        argv.reserve(arguments.size() + 1U);
        for (auto &argument : arguments) {
            argv.push_back(const_cast<char *>(argument.c_str()));
        }
        argv.push_back(nullptr);

        const auto start = std::chrono::steady_clock::now();
        pid_t child = 0;
        const int spawn_error =
            ::posix_spawnp(&child, arguments.front().c_str(),
                           file_actions.get(), attributes.get(), argv.data(),
                           environ);
        if (spawn_error != 0) {
            throw std::runtime_error(
                "could not launch child: " +
                std::string(std::strerror(spawn_error)));
        }

        ChildGuard child_guard(child);
        const WaitResult wait_result =
            waitForChild(child, config.timeout);
        child_guard.release();
        const auto end = std::chrono::steady_clock::now();

        if (wait_result.timed_out) {
            throw std::runtime_error(
                "timed out after " +
                std::to_string(config.timeout.count()) +
                " ms; child terminated with " +
                describeStatus(wait_result.status));
        }
        if (!WIFEXITED(wait_result.status) ||
            WEXITSTATUS(wait_result.status) != 0) {
            throw std::runtime_error(
                "child terminated with " +
                describeStatus(wait_result.status));
        }
        if (!pathExists(artefacts.resultPath())) {
            throw std::runtime_error(
                "child exited successfully but did not create a result");
        }
        if (!pathExists(artefacts.metricsPath())) {
            throw std::runtime_error(
                "child exited successfully but did not create metrics");
        }

        InferenceExecution execution;
        execution.elapsed_ms =
            std::chrono::duration<double, std::milli>(end - start).count();
        return execution;
    } catch (const std::exception &error) {
        processError(config, artefacts, error.what());
    }
}

std::string InferenceProcess::diagnostics(
    const InferenceArtefacts &artefacts, const std::size_t max_bytes) {
    std::ostringstream output;
    output << "--- child stdout (bounded) ---\n"
           << readTail(artefacts.stdoutPath(), max_bytes)
           << "\n--- child stderr (bounded) ---\n"
           << readTail(artefacts.stderrPath(), max_bytes)
           << "\n--- end child diagnostics ---";
    return output.str();
}

} // namespace image
