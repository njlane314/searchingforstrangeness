#ifndef INFERENCEPRODUCTION_H
#define INFERENCEPRODUCTION_H

#include "Products/ImageFeatures.h"

#include <chrono>
#include <cerrno>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <fcntl.h>
#include <spawn.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

extern char **environ;

namespace image {

class InferenceProduction {
public:
    struct Perf {
        double t_write_req_ms{0.0};
        double t_exec_total_ms{0.0};
        double t_read_resp_ms{0.0};
        double t_child_total_ms{0.0};
        double t_child_setup_ms{0.0};
        double t_child_infer_ms{0.0};
        double t_child_post_ms{0.0};
        double child_max_rss_mb{0.0};
    };

    struct Result {
        std::vector<float> cls;
        Perf perf;
    };

    static Result runInference(const std::vector<ImageFeatures> &detector_images,
                               const std::string &absolute_scratch_dir,
                               const std::string &arch,
                               const std::string &weights_file,
                               const std::string &inference_wrapper,
                               const std::string &runtime_base_dir);
};

inline std::string joinPath(std::string a, const std::string &b) {
    if (a.empty()) return b;
    if (!a.empty() && a.back() != '/') a.push_back('/');
    return a + b;
}

inline std::string makeAbsolutePath(const std::string &path) {
    if (path.empty())
        return path;

    if (path.front() == '/' ||
        (path.size() > 1 && path[1] == ':') ||
        path.find("://") != std::string::npos)
        return path;

    char resolved[4096];
    if (realpath(path.c_str(), resolved))
        return std::string(resolved);

    char cwd[4096];
    if (getcwd(cwd, sizeof(cwd)))
        return joinPath(std::string(cwd), path);

    return path;
}

inline bool isLocalFilesystemPath(const std::string &path) {
    return !path.empty() &&
           path.front() == '/' &&
           path.find("://") == std::string::npos;
}

inline std::string parentDir(const std::string &path) {
    if (!isLocalFilesystemPath(path))
        return {};

    auto pos = path.rfind('/');
    if (pos == std::string::npos)
        return {};
    if (pos == 0)
        return "/";
    return path.substr(0, pos);
}

inline void appendUniqueBind(std::vector<std::string> &binds,
                             std::set<std::string> &seen,
                             const std::string &path) {
    if (path.empty())
        return;
    if (seen.insert(path).second)
        binds.push_back(path);
}

inline std::string joinArgsForLog(const std::vector<std::string> &args) {
    std::ostringstream os;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i)
            os << ' ';
        os << args[i];
    }
    return os.str();
}

namespace _binary_io {
    struct RequestHeader {
        char magic[4];
        uint32_t version;
        uint32_t C;
        uint32_t H;
        uint32_t W;
        uint32_t F;
    };

    struct PlaneHeader {
        uint64_t nnz;
    };

    struct ResultHeader {
        char magic[4];
        uint32_t version;
        uint32_t K;
        uint32_t reserved;
        uint64_t cls_bytes;
    };

    static_assert(sizeof(RequestHeader) == 24, "Unexpected RequestHeader size");
    static_assert(sizeof(PlaneHeader) == 8, "Unexpected PlaneHeader size");
    static_assert(sizeof(ResultHeader) == 24, "Unexpected ResultHeader size");

    inline void write_sparse_plane(std::ofstream &ofs,
                                   const ImageFeatures &plane) {
        if (plane.feature_dim == 0) {
            throw std::runtime_error("Sparse plane feature_dim must be non-zero.");
        }
        if ((plane.coords.size() % 2) != 0) {
            throw std::runtime_error("Sparse plane coords payload must contain row/col pairs.");
        }
        std::size_t const nnz = plane.coords.size() / 2;
        if (plane.features.size() != nnz * static_cast<std::size_t>(plane.feature_dim)) {
            throw std::runtime_error("Sparse plane payload is inconsistent.");
        }

        PlaneHeader h{static_cast<uint64_t>(nnz)};
        ofs.write(reinterpret_cast<const char *>(&h), sizeof(h));
        if (h.nnz) {
            ofs.write(reinterpret_cast<const char *>(plane.coords.data()),
                      sizeof(int32_t) * plane.coords.size());
            ofs.write(reinterpret_cast<const char *>(plane.features.data()),
                      sizeof(float) * plane.features.size());
        }
        if (!ofs)
            throw std::runtime_error("Short write while writing sparse plane payload.");
    }

    inline void write_sparse_request(const std::string &path,
                                     const ImageFeatures &u,
                                     const ImageFeatures &v,
                                     const ImageFeatures &w) {
        if (u.width != v.width || u.width != w.width ||
            u.height != v.height || u.height != w.height) {
            throw std::runtime_error("Sparse inference request requires matching plane dimensions.");
        }
        if (u.feature_dim != v.feature_dim || u.feature_dim != w.feature_dim) {
            throw std::runtime_error("Sparse inference request requires matching plane feature dimensions.");
        }

        std::ofstream ofs(path, std::ios::binary | std::ios::trunc);
        if (!ofs)
            throw std::runtime_error("Cannot open " + path);

        RequestHeader header{{'I', 'A', 'S', 'P'}, 2, 3, u.height, u.width, u.feature_dim};
        ofs.write(reinterpret_cast<const char *>(&header), sizeof(header));
        write_sparse_plane(ofs, u);
        write_sparse_plane(ofs, v);
        write_sparse_plane(ofs, w);
        if (!ofs)
            throw std::runtime_error("Short write to " + path);
    }
}

inline InferenceProduction::Result InferenceProduction::runInference(
    const std::vector<ImageFeatures> &detector_images,
    const std::string &absolute_scratch_dir,
    const std::string &arch, const std::string &weights_file,
    const std::string &inference_wrapper,
    const std::string &runtime_base_dir) {
    using std::string;

    if (detector_images.size() < 3) {
        throw art::Exception(art::errors::LogicError)
            << "Need at least 3 detector images (U,V,W)";
    }

    auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    int pid  = ::getpid();
    std::ostringstream base;
    base << "ia_" << pid << "_" << now;

    string req_bin       = joinPath(absolute_scratch_dir, base.str() + "_planes.bin");
    string result_bin    = joinPath(absolute_scratch_dir, base.str() + "_results.bin");
    string meta_txt      = result_bin + ".meta";
    string script_stdout = joinPath(absolute_scratch_dir, base.str() + "_py_script.out");
    string script_stderr = joinPath(absolute_scratch_dir, base.str() + "_py_script.err");

    auto t0 = std::chrono::steady_clock::now();
    const auto &U = detector_images[0];
    const auto &V = detector_images[1];
    const auto &W = detector_images[2];
    _binary_io::write_sparse_request(req_bin, U, V, W);
    auto t1 = std::chrono::steady_clock::now();

    string container = "/cvmfs/uboone.opensciencegrid.org/containers/lantern_v2_me_06_03_prod";

    std::string runtime_dir = runtime_base_dir;
    if (runtime_dir.empty()) {
        auto pos = inference_wrapper.rfind("/scripts/");
        runtime_dir =
            (pos == std::string::npos) ? std::string{} : inference_wrapper.substr(0, pos);
    }

    std::string abs_runtime_dir = makeAbsolutePath(runtime_dir);
    std::string abs_inference_wrapper = makeAbsolutePath(inference_wrapper);
    std::string abs_weights_file = makeAbsolutePath(weights_file);

    std::vector<std::string> binds;
    std::set<std::string> bind_seen;
    appendUniqueBind(binds, bind_seen, absolute_scratch_dir);
    appendUniqueBind(binds, bind_seen, abs_runtime_dir);
    appendUniqueBind(binds, bind_seen, parentDir(abs_inference_wrapper));
    appendUniqueBind(binds, bind_seen, parentDir(abs_weights_file));

    bool need_cvmfs = true;
    if (const char *bind_env = std::getenv("APPTAINER_BINDPATH")) {
        if (std::string(bind_env).find("/cvmfs") != std::string::npos)
            need_cvmfs = false;
    } else {
        std::ifstream mounts("/proc/mounts");
        std::string line;
        while (std::getline(mounts, line)) {
            if (line.find(" /cvmfs ") != std::string::npos) {
                need_cvmfs = false;
                break;
            }
        }
    }
    if (need_cvmfs)
        binds.insert(binds.begin(), "/cvmfs");

    std::ostringstream bind_csv;
    for (size_t i = 0; i < binds.size(); ++i) {
        if (i)
            bind_csv << ",";
        bind_csv << binds[i];
    }

    std::vector<std::string> args{
        "apptainer",
        "exec",
        "--cleanenv",
        "--bind",
        bind_csv.str(),
        container,
        "/bin/bash",
        abs_inference_wrapper,
        "--in",
        req_bin,
        "--out",
        result_bin,
        "--metrics",
        meta_txt,
        "--W",
        std::to_string(detector_images[0].width),
        "--H",
        std::to_string(detector_images[0].height),
        "--arch",
        arch,
        "--weights",
        abs_weights_file,
    };

    mf::LogInfo("InferenceProduction")
        << "Executing inference: " << joinArgsForLog(args);

    auto start = std::chrono::steady_clock::now();
    posix_spawn_file_actions_t file_actions;
    if (int rc = posix_spawn_file_actions_init(&file_actions); rc != 0) {
        throw art::Exception(art::errors::LogicError)
            << "Failed to initialize spawn file actions: " << std::strerror(rc);
    }

    bool file_actions_ready = true;
    auto destroy_file_actions = [&]() {
        if (file_actions_ready) {
            posix_spawn_file_actions_destroy(&file_actions);
            file_actions_ready = false;
        }
    };

    auto add_open = [&](int fd, const std::string &path) {
        int rc = posix_spawn_file_actions_addopen(
            &file_actions, fd, path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (rc != 0) {
            destroy_file_actions();
            throw art::Exception(art::errors::LogicError)
                << "Failed to redirect fd " << fd << " to " << path
                << ": " << std::strerror(rc);
        }
    };
    add_open(STDOUT_FILENO, script_stdout);
    add_open(STDERR_FILENO, script_stderr);

    std::vector<char *> argv;
    argv.reserve(args.size() + 1);
    for (auto &arg : args)
        argv.push_back(const_cast<char *>(arg.c_str()));
    argv.push_back(nullptr);

    pid_t child_pid = 0;
    int spawn_rc = posix_spawnp(&child_pid, args.front().c_str(),
                                &file_actions, nullptr, argv.data(), environ);
    destroy_file_actions();
    if (spawn_rc != 0) {
        throw art::Exception(art::errors::LogicError)
            << "Failed to launch inference process: " << std::strerror(spawn_rc);
    }

    int status = 0;
    if (waitpid(child_pid, &status, 0) < 0) {
        throw art::Exception(art::errors::LogicError)
            << "waitpid failed for inference process: " << std::strerror(errno);
    }
    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();

    struct stat temp_stat;
    bool exited_cleanly = WIFEXITED(status) && (WEXITSTATUS(status) == 0);
    bool success = exited_cleanly && (stat(result_bin.c_str(), &temp_stat) == 0);
    if (!success) {
        std::ifstream error_stream(script_stderr);
        std::string error_message(
            (std::istreambuf_iterator<char>(error_stream)), {});
        std::ostringstream status_message;
        if (WIFEXITED(status)) {
            status_message << "exit code " << WEXITSTATUS(status);
        } else if (WIFSIGNALED(status)) {
            status_message << "signal " << WTERMSIG(status);
        } else {
            status_message << "status " << status;
        }
        throw art::Exception(art::errors::LogicError)
            << "Inference script failed with " << status_message.str()
            << "\n--- Script stderr ---\n"
            << error_message << "\n--- End Script stderr ---";
    }

    auto t2 = std::chrono::steady_clock::now();
    Result out;
    std::ifstream ifs(result_bin, std::ios::binary);
    if (!ifs) {
        throw art::Exception(art::errors::LogicError)
            << "Could not open result file: " << result_bin;
    }
    _binary_io::ResultHeader h{};
    ifs.read(reinterpret_cast<char *>(&h), sizeof(h));
    if (!ifs || std::string(h.magic, h.magic + 4) != "IAOK" || h.version != 1) {
        throw art::Exception(art::errors::LogicError)
            << "Bad result header in " << result_bin;
    }
    if (h.K && h.cls_bytes == h.K * sizeof(float)) {
        out.cls.resize(h.K);
        ifs.read(reinterpret_cast<char *>(out.cls.data()), h.cls_bytes);
    }

    auto t3 = std::chrono::steady_clock::now();

    out.perf.t_write_req_ms  = std::chrono::duration<double, std::milli>(t1 - t0).count();
    out.perf.t_exec_total_ms = duration * 1000.0;
    out.perf.t_read_resp_ms  = std::chrono::duration<double, std::milli>(t3 - t2).count();

    {
        std::ifstream m(meta_txt);
        if (m) {
            std::string line;
            while (std::getline(m, line)) {
                auto eq = line.find('=');
                if (eq == std::string::npos) continue;
                std::string key = line.substr(0, eq);
                std::string val = line.substr(eq + 1);

                char *endp = nullptr;
                double d = std::strtod(val.c_str(), &endp);
                bool numeric = (endp != val.c_str());

                if (key == "t_total_ms" && numeric)
                    out.perf.t_child_total_ms = d;
                else if (key == "t_setup_ms" && numeric)
                    out.perf.t_child_setup_ms = d;
                else if (key == "t_infer_ms" && numeric)
                    out.perf.t_child_infer_ms = d;
                else if (key == "t_post_ms" && numeric)
                    out.perf.t_child_post_ms = d;
                else if (key == "max_rss_mb" && numeric)
                    out.perf.child_max_rss_mb = d;
            }
        }
    }

    mf::LogInfo("InferenceProduction")
        << "Inference time: " << duration << " seconds";
    if (!out.cls.empty())
        mf::LogInfo("InferenceProduction") << "First class score: " << out.cls.front();
    else
        mf::LogWarning("InferenceProduction")
            << "No class scores returned; out.cls is empty";
    mf::LogInfo("InferenceProduction")
        << "t_write_req_ms=" << out.perf.t_write_req_ms
        << " t_exec_total_ms=" << out.perf.t_exec_total_ms
        << " t_read_resp_ms=" << out.perf.t_read_resp_ms
        << " child(total/setup/infer/post)_ms="
        << out.perf.t_child_total_ms << "/"
        << out.perf.t_child_setup_ms << "/"
        << out.perf.t_child_infer_ms << "/"
        << out.perf.t_child_post_ms;

    std::remove(req_bin.c_str());
    std::remove(result_bin.c_str());
    std::remove(script_stdout.c_str());
    std::remove(script_stderr.c_str());
    std::remove(meta_txt.c_str());
    return out;
}

}

#endif
