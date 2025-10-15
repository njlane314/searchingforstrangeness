#ifndef INFERENCEENGINE_H
#define INFERENCEENGINE_H

#include "Imaging/Image.h"

#include <chrono>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
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
class InferenceEngine {
  public:
    struct Result {
      std::vector<float> cls;
      uint32_t segW{0}, segH{0};
      std::vector<uint8_t> seg_u, seg_v, seg_w;
      std::vector<float> conf_u, conf_v, conf_w;
    };

    static Result runInferenceDetailed(const std::vector<Image<float>> &detector_images,
                                       const std::string &absolute_scratch_dir,
                                       const std::string &work_dir,
                                       const std::string &arch,
                                       const std::string &weights_file,
                                       const std::string &inference_wrapper,
                                       const std::string &assets_base_dir,
                                       bool want_cls = true,
                                       bool want_seg = true);

    static float runInference(const std::vector<Image<float>> &detector_images,
                              const std::string &absolute_scratch_dir,
                              const std::string &work_dir,
                              const std::string &arch,
                              const std::string &weights_file,
                              const std::string &inference_wrapper,
                              const std::string &assets_base_dir) {
      auto r = runInferenceDetailed(detector_images, absolute_scratch_dir, work_dir,
                                    arch, weights_file, inference_wrapper, assets_base_dir,
                                    true, false);
      return r.cls.empty() ? std::numeric_limits<float>::quiet_NaN() : r.cls.front();
    }
};

inline std::string joinPath(std::string a, const std::string &b) {
    if (a.empty()) return b;
    if (!a.empty() && a.back() != '/') a.push_back('/');
    return a + b;
}

namespace _binary_io {
  struct ResultHeader {
    char     magic[4];
    uint32_t version;
    uint32_t K;
    uint32_t segW;
    uint32_t segH;
    uint32_t has_conf;
    uint64_t cls_bytes;
    uint64_t seg_bytes;
    uint64_t conf_bytes;
  };

  inline void write_chw_f32(const std::string& path,
                            const std::vector<float>& u,
                            const std::vector<float>& v,
                            const std::vector<float>& w) {
    std::ofstream ofs(path, std::ios::binary | std::ios::trunc);
    if (!ofs) throw std::runtime_error("Cannot open " + path);
    ofs.write(reinterpret_cast<const char*>(u.data()), sizeof(float)*u.size());
    ofs.write(reinterpret_cast<const char*>(v.data()), sizeof(float)*v.size());
    ofs.write(reinterpret_cast<const char*>(w.data()), sizeof(float)*w.size());
    if (!ofs) throw std::runtime_error("Short write to " + path);
  }
}

inline InferenceEngine::Result InferenceEngine::runInferenceDetailed(
    const std::vector<Image<float>> &detector_images,
    const std::string &absolute_scratch_dir, const std::string &work_dir,
    const std::string &arch, const std::string &weights_file,
    const std::string &inference_wrapper,
    const std::string &assets_base_dir,
    bool want_cls, bool want_seg) {
    using std::string;
    (void)work_dir;

    auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    int pid  = ::getpid();
    std::ostringstream base;
    base << "ia_" << pid << "_" << now;

    string req_bin       = joinPath(absolute_scratch_dir, base.str() + "_planes.bin");
    string result_bin    = joinPath(absolute_scratch_dir, base.str() + "_results.bin");
    string script_stdout = joinPath(absolute_scratch_dir, base.str() + "_py_script.out");
    string script_stderr = joinPath(absolute_scratch_dir, base.str() + "_py_script.err");

    const auto& U = detector_images[0].data();
    const auto& V = detector_images[1].data();
    const auto& W = detector_images[2].data();
    _binary_io::write_chw_f32(req_bin, U, V, W);

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
        << "--W " << detector_images[0].width() << " "
        << "--H " << detector_images[0].height() << " "
        << "--arch " << arch << " "
        << "--weights " << weights_file << " "
        << "--want-cls " << (want_cls ? 1 : 0) << " "
        << "--want-seg " << (want_seg ? 1 : 0)
        << " > " << script_stdout << " 2> " << script_stderr;

    mf::LogInfo("InferenceEngine") << "Executing inference: " << cmd.str();

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

    Result out;
    std::ifstream ifs(result_bin, std::ios::binary);
    if (!ifs) {
      throw art::Exception(art::errors::LogicError)
        << "Could not open result file: " << result_bin;
    }
    _binary_io::ResultHeader h{};
    ifs.read(reinterpret_cast<char*>(&h), sizeof(h));
    if (!ifs || std::string(h.magic, h.magic+4) != "IAOK" || h.version != 1) {
      throw art::Exception(art::errors::LogicError)
        << "Bad result header in " << result_bin;
    }
    if (h.K && h.cls_bytes == h.K * sizeof(float)) {
      out.cls.resize(h.K);
      ifs.read(reinterpret_cast<char*>(out.cls.data()), h.cls_bytes);
    }
    if (h.segW && h.segH) {
      out.segW = h.segW; out.segH = h.segH;
      const size_t npix = size_t(out.segW) * out.segH;
      if (h.seg_bytes != 3 * npix) {
        throw art::Exception(art::errors::LogicError)
          << "seg_bytes mismatch in " << result_bin;
      }
      out.seg_u.resize(npix); out.seg_v.resize(npix); out.seg_w.resize(npix);
      ifs.read(reinterpret_cast<char*>(out.seg_u.data()), npix);
      ifs.read(reinterpret_cast<char*>(out.seg_v.data()), npix);
      ifs.read(reinterpret_cast<char*>(out.seg_w.data()), npix);
      if (h.has_conf) {
        if (h.conf_bytes != 3 * npix * sizeof(float)) {
          throw art::Exception(art::errors::LogicError)
            << "conf_bytes mismatch in " << result_bin;
        }
        out.conf_u.resize(npix); out.conf_v.resize(npix); out.conf_w.resize(npix);
        ifs.read(reinterpret_cast<char*>(out.conf_u.data()), npix*sizeof(float));
        ifs.read(reinterpret_cast<char*>(out.conf_v.data()), npix*sizeof(float));
        ifs.read(reinterpret_cast<char*>(out.conf_w.data()), npix*sizeof(float));
      }
    }

    mf::LogInfo("InferenceEngine")
        << "Inference time: " << duration << " seconds";
    if (!out.cls.empty())
      mf::LogInfo("InferenceEngine") << "First class score: " << out.cls.front();

    std::remove(req_bin.c_str());
    std::remove(result_bin.c_str());
    std::remove(script_stdout.c_str());
    std::remove(script_stderr.c_str());

    return out;
}

} // namespace image

#endif
