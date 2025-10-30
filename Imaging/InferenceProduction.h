#ifndef INFERENCEPRODUCTION_H
#define INFERENCEPRODUCTION_H

#include "Products/ImageProducts.h"

#include <chrono>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <cetlib_except/exception.h>
#include <canvas/Utilities/Exception.h>
#include <messagefacility/MessageLogger/MessageLogger.h>

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

        static Result runInference(const std::vector<PlaneImage> &detector_images,
                                   const std::string &absolute_scratch_dir,
                                   const std::string &work_dir,
                                   const std::string &arch,
                                   const std::string &weights_file,
                                   const std::string &inference_wrapper,
                                   const std::string &assets_base_dir);
    };

inline std::string joinPath(std::string a, const std::string &b) {
    if (a.empty()) return b;
    if (!a.empty() && a.back() != '/') a.push_back('/');
    return a + b;
}

namespace _binary_io {
    struct ResultHeader {
        char magic[4];
        uint32_t version;
        uint32_t K;
        uint32_t reserved;
        uint64_t cls_bytes;
    };

    static_assert(sizeof(ResultHeader) == 24, "Unexpected ResultHeader size");

    inline void write_chw_f32(const std::string &path,
                              const std::vector<float> &u,
                              const std::vector<float> &v,
                              const std::vector<float> &w) {
        std::ofstream ofs(path, std::ios::binary | std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open " + path);
        ofs.write(reinterpret_cast<const char *>(u.data()), sizeof(float) * u.size());
        ofs.write(reinterpret_cast<const char *>(v.data()), sizeof(float) * v.size());
        ofs.write(reinterpret_cast<const char *>(w.data()), sizeof(float) * w.size());
        if (!ofs) throw std::runtime_error("Short write to " + path);
    }
}

inline InferenceProduction::Result InferenceProduction::runInference(
    const std::vector<PlaneImage> &detector_images,
    const std::string &absolute_scratch_dir, const std::string &work_dir,
    const std::string &arch, const std::string &weights_file,
    const std::string &inference_wrapper,
    const std::string &assets_base_dir) {
    using std::string;
    (void)work_dir;

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

    // NOTE: wrapper will also write a flat float32 sidecar "<result_bin>.feat.f32"

    auto t0 = std::chrono::steady_clock::now();
    const auto &U = detector_images[0].adc;
    const auto &V = detector_images[1].adc;
    const auto &W = detector_images[2].adc;
    _binary_io::write_chw_f32(req_bin, U, V, W);
    auto t1 = std::chrono::steady_clock::now();

    string container = "/cvmfs/uboone.opensciencegrid.org/containers/lantern_v2_me_06_03_prod";

    std::string assets_dir = assets_base_dir;
    if (assets_dir.empty()) {
        auto pos = inference_wrapper.rfind("/scripts/");
        assets_dir = (pos == std::string::npos) ? std::string{} : inference_wrapper.substr(0, pos);
    }

    std::vector<std::string> binds;
    if (!assets_dir.empty())
        binds.push_back(assets_dir);
    binds.push_back(absolute_scratch_dir);

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

    std::ostringstream cmd;
    cmd << "apptainer exec --cleanenv --bind " << bind_csv.str() << " "
        << container << " "
        << "/bin/bash " << inference_wrapper << " "
        << "--in " << req_bin << " "
        << "--out " << result_bin << " "
        << "--metrics " << meta_txt << " "
        << "--W " << detector_images[0].width << " "
        << "--H " << detector_images[0].height << " "
        << "--arch " << arch << " "
        << "--weights " << weights_file << " "
        << " > " << script_stdout << " 2> " << script_stderr;

    mf::LogInfo("InferenceProduction") << "Executing inference: " << cmd.str();

    auto start = std::chrono::steady_clock::now();
    int code = std::system(cmd.str().c_str());
    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();

    struct stat temp_stat;
    bool success = (code == 0) && (stat(result_bin.c_str(), &temp_stat) == 0);
    if (!success) {
        std::ifstream error_stream(script_stderr);
        std::string error_message(
            (std::istreambuf_iterator<char>(error_stream)), {});
        throw art::Exception(art::errors::LogicError)
            << "Inference script failed with exit code " << code
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

    // Parse metrics text
    {
        std::ifstream m(meta_txt);
        if (m) {
            std::string line;
            while (std::getline(m, line)) {
                auto eq = line.find('=');
                if (eq == std::string::npos) continue;
                std::string key = line.substr(0, eq);
                std::string val = line.substr(eq + 1);

                // numeric parse
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
    // keep the .feat.f32 if you want to inspect offline (we do not delete here)

    return out;
}

}

#endif
