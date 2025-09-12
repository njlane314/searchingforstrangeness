#ifndef INFERENCEENGINE_H
#define INFERENCEENGINE_H

#include "Imaging/Image.h"
#include "Common/NpyUtils.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <cstdio>
#include <cetlib_except/exception.h>
#include <canvas/Utilities/Exception.h>
#include <messagefacility/MessageLogger/MessageLogger.h>

namespace analysis {
class InferenceEngine {
  public:
    static float runInference(const std::vector<Image<float>> &detector_images,
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

inline float InferenceEngine::runInference(
    const std::vector<Image<float>> &detector_images,
    const std::string &absolute_scratch_dir, const std::string &work_dir,
    const std::string &arch, const std::string &weights_file,
    const std::string &inference_wrapper,
    const std::string &assets_base_dir) {
    using std::string;
    (void)work_dir;
    string npy_in = joinPath(absolute_scratch_dir, "detector_images.npy");
    string temp_out = joinPath(absolute_scratch_dir, "temp_test_out.txt");
    string script_stdout = joinPath(absolute_scratch_dir, "py_script.out");
    string script_stderr = joinPath(absolute_scratch_dir, "py_script.err");

    std::vector<std::vector<float>> images = {
        detector_images[0].data(),
        detector_images[1].data(),
        detector_images[2].data()};
    common::save_npy_f32_2d(npy_in, images);

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
        << "--npy " << npy_in << " "
        << "--output " << temp_out << " "
        << "--arch " << arch << " "
        << "--weights " << weights_file
        << " > " << script_stdout << " 2> " << script_stderr;

    mf::LogInfo("InferenceEngine") << "Executing inference: " << cmd.str();

    auto start = std::chrono::steady_clock::now();
    int code = std::system(cmd.str().c_str());
    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();

    struct stat temp_stat;
    bool success = (code == 0) && (stat(temp_out.c_str(), &temp_stat) == 0);
    if (!success) {
        std::ifstream error_stream(script_stderr);
        std::string error_message(
            (std::istreambuf_iterator<char>(error_stream)), {});
        throw art::Exception(art::errors::LogicError)
            << "Inference script failed with exit code " << code
            << "\n--- Script stderr ---\n"
            << error_message << "\n--- End Script stderr ---";
    }

    std::ifstream result_stream(temp_out.c_str());
    if (!result_stream) {
        throw art::Exception(art::errors::LogicError)
            << "Could not open temporary result file: " << temp_out;
    }
    float score;
    result_stream >> score;

    mf::LogInfo("InferenceEngine")
        << "Inference time: " << duration << " seconds";
    mf::LogInfo("InferenceEngine") << "Predicted score: " << score;

    std::remove(npy_in.c_str());
    std::remove(temp_out.c_str());
    std::remove(script_stdout.c_str());
    std::remove(script_stderr.c_str());

    return score;
}

} // namespace analysis

#endif // INFERENCEENGINE_H
