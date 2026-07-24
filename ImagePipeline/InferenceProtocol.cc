#include "ImagePipeline/InferenceProtocol.h"

#include <array>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace image {
namespace {

static_assert(sizeof(float) == sizeof(uint32_t),
              "IASP/IAOK require 32-bit IEEE-754 floats");
static_assert(std::numeric_limits<float>::is_iec559,
              "IASP/IAOK require IEEE-754 floats");

[[noreturn]] void protocolError(const std::string &path,
                                const std::string &message) {
    throw std::runtime_error("Inference protocol error in " + path + ": " +
                             message);
}

void writeExact(std::ofstream &stream, const char *data, std::size_t bytes,
                const std::string &path) {
    stream.write(data, static_cast<std::streamsize>(bytes));
    if (!stream) {
        protocolError(path, "short write");
    }
}

void writeU32(std::ofstream &stream, uint32_t value,
              const std::string &path) {
    const std::array<char, 4> bytes{{
        static_cast<char>(value & 0xffU),
        static_cast<char>((value >> 8U) & 0xffU),
        static_cast<char>((value >> 16U) & 0xffU),
        static_cast<char>((value >> 24U) & 0xffU),
    }};
    writeExact(stream, bytes.data(), bytes.size(), path);
}

void writeU64(std::ofstream &stream, uint64_t value,
              const std::string &path) {
    std::array<char, 8> bytes{};
    for (std::size_t index = 0; index < bytes.size(); ++index) {
        bytes[index] =
            static_cast<char>((value >> (8U * index)) & uint64_t{0xffU});
    }
    writeExact(stream, bytes.data(), bytes.size(), path);
}

void writeI32(std::ofstream &stream, int32_t value,
              const std::string &path) {
    uint32_t bits = 0U;
    static_assert(sizeof(bits) == sizeof(value), "Unexpected int32_t size");
    std::memcpy(&bits, &value, sizeof(bits));
    writeU32(stream, bits, path);
}

void writeFloat(std::ofstream &stream, float value,
                const std::string &path) {
    if (!std::isfinite(value)) {
        protocolError(path, "cannot serialise a non-finite feature");
    }
    uint32_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    writeU32(stream, bits, path);
}

void writePlane(std::ofstream &stream, const ImageFeatures &plane,
                const std::string &path) {
    const std::size_t active_pixels = plane.coords.size() / 2U;
    if (active_pixels >
        static_cast<std::size_t>(std::numeric_limits<uint64_t>::max())) {
        protocolError(path, "plane active-pixel count exceeds uint64");
    }

    writeU64(stream, static_cast<uint64_t>(active_pixels), path);
    for (const auto coordinate : plane.coords) {
        writeI32(stream, coordinate, path);
    }
    for (const auto feature : plane.features) {
        writeFloat(stream, feature, path);
    }
}

void readExact(std::ifstream &stream, char *data, std::size_t bytes,
               const std::string &path, const std::string &field) {
    stream.read(data, static_cast<std::streamsize>(bytes));
    if (stream.gcount() != static_cast<std::streamsize>(bytes)) {
        protocolError(path, "truncated " + field);
    }
}

uint32_t readU32(std::ifstream &stream, const std::string &path,
                 const std::string &field) {
    std::array<unsigned char, 4> bytes{};
    readExact(stream, reinterpret_cast<char *>(bytes.data()), bytes.size(),
              path, field);
    return static_cast<uint32_t>(bytes[0]) |
           (static_cast<uint32_t>(bytes[1]) << 8U) |
           (static_cast<uint32_t>(bytes[2]) << 16U) |
           (static_cast<uint32_t>(bytes[3]) << 24U);
}

uint64_t readU64(std::ifstream &stream, const std::string &path,
                 const std::string &field) {
    std::array<unsigned char, 8> bytes{};
    readExact(stream, reinterpret_cast<char *>(bytes.data()), bytes.size(),
              path, field);

    uint64_t value = 0U;
    for (std::size_t index = 0; index < bytes.size(); ++index) {
        value |= static_cast<uint64_t>(bytes[index]) << (8U * index);
    }
    return value;
}

float readFloat(std::ifstream &stream, const std::string &path,
                const std::string &field) {
    const uint32_t bits = readU32(stream, path, field);
    float value = 0.f;
    std::memcpy(&value, &bits, sizeof(value));
    if (!std::isfinite(value)) {
        protocolError(path, field + " is non-finite");
    }
    return value;
}

double parseMetricValue(const std::string &path, const std::string &key,
                        const std::string &text) {
    if (text.empty() ||
        std::isspace(static_cast<unsigned char>(text.front())) ||
        std::isspace(static_cast<unsigned char>(text.back()))) {
        protocolError(path, "invalid value for metrics key '" + key + "'");
    }

    errno = 0;
    char *end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
        !std::isfinite(value)) {
        protocolError(path, "invalid numeric value for metrics key '" + key +
                                "': " + text);
    }
    return value;
}

} // namespace

void InferenceProtocol::writeRequest(
    const std::string &path, const ValidatedPlaneTriplet &planes) {
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream) {
        protocolError(path, "cannot open request for writing");
    }

    constexpr std::array<char, 4> magic{{'I', 'A', 'S', 'P'}};
    writeExact(stream, magic.data(), magic.size(), path);
    writeU32(stream, requestVersion(), path);
    writeU32(stream, 3U, path);
    writeU32(stream, planes.height(), path);
    writeU32(stream, planes.width(), path);
    writeU32(stream, planes.featureDimension(), path);

    // Protocol plane order is a stable part of the model contract.
    writePlane(stream, planes.u(), path);
    writePlane(stream, planes.v(), path);
    writePlane(stream, planes.w(), path);

    stream.close();
    if (!stream) {
        protocolError(path, "failed to finalise request");
    }
}

std::vector<float>
InferenceProtocol::readResult(const std::string &path,
                              const std::size_t max_logits) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) {
        protocolError(path, "cannot open result");
    }

    std::array<char, 4> magic{};
    readExact(stream, magic.data(), magic.size(), path, "magic");
    if (magic != std::array<char, 4>{{'I', 'A', 'O', 'K'}}) {
        protocolError(path, "bad result magic");
    }

    const uint32_t version = readU32(stream, path, "version");
    if (version != resultVersion()) {
        protocolError(path, "unsupported IAOK version " +
                                std::to_string(version));
    }

    const uint32_t count = readU32(stream, path, "logit count");
    const uint32_t reserved = readU32(stream, path, "reserved field");
    const uint64_t payload_bytes =
        readU64(stream, path, "logit payload byte count");

    if (reserved != 0U) {
        protocolError(path, "reserved result field must be zero");
    }
    if (static_cast<std::size_t>(count) > max_logits) {
        protocolError(path, "result logit count " + std::to_string(count) +
                                " exceeds configured cap " +
                                std::to_string(max_logits));
    }

    const uint64_t expected_bytes =
        static_cast<uint64_t>(count) * sizeof(float);
    if (payload_bytes != expected_bytes) {
        protocolError(path,
                      "declared payload byte count " +
                          std::to_string(payload_bytes) +
                          " does not match K*sizeof(float) (" +
                          std::to_string(expected_bytes) + ")");
    }

    std::vector<float> logits;
    logits.reserve(count);
    for (uint32_t index = 0; index < count; ++index) {
        logits.push_back(
            readFloat(stream, path, "logit " + std::to_string(index)));
    }

    char trailing = '\0';
    stream.read(&trailing, 1);
    if (stream.gcount() != 0) {
        protocolError(path, "unexpected trailing bytes after result payload");
    }
    if (!stream.eof()) {
        protocolError(path, "failed while checking result file boundary");
    }

    return logits;
}

InferenceChildMetrics
InferenceProtocol::readMetrics(const std::string &path) {
    std::ifstream stream(path);
    if (!stream) {
        protocolError(path, "cannot open metrics");
    }

    std::map<std::string, double> values;
    std::string line;
    std::size_t line_number = 0U;
    while (std::getline(stream, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        const auto separator = line.find('=');
        if (separator == std::string::npos ||
            separator == 0U ||
            separator + 1U >= line.size() ||
            line.find('=', separator + 1U) != std::string::npos) {
            protocolError(path, "malformed metrics line " +
                                    std::to_string(line_number));
        }

        const std::string key = line.substr(0U, separator);
        const std::string text = line.substr(separator + 1U);
        if (key != "t_total_ms" && key != "t_setup_ms" &&
            key != "t_infer_ms" && key != "t_post_ms" &&
            key != "max_rss_mb" && key != "seed" && key != "thr") {
            protocolError(path, "unknown metrics key '" + key + "'");
        }
        if (values.count(key) != 0U) {
            protocolError(path, "duplicate metrics key '" + key + "'");
        }
        values.emplace(key, parseMetricValue(path, key, text));
    }
    if (stream.bad()) {
        protocolError(path, "I/O failure while reading metrics");
    }

    constexpr std::array<const char *, 5> required{{
        "t_total_ms",
        "t_setup_ms",
        "t_infer_ms",
        "t_post_ms",
        "max_rss_mb",
    }};
    for (const auto *key : required) {
        const auto found = values.find(key);
        if (found == values.end()) {
            protocolError(path, "missing required metrics key '" +
                                    std::string(key) + "'");
        }
        if (found->second < 0.0) {
            protocolError(path, "metrics key '" + std::string(key) +
                                    "' must be non-negative");
        }
    }

    InferenceChildMetrics metrics;
    metrics.total_ms = values.at("t_total_ms");
    metrics.setup_ms = values.at("t_setup_ms");
    metrics.infer_ms = values.at("t_infer_ms");
    metrics.post_ms = values.at("t_post_ms");
    metrics.max_rss_mb = values.at("max_rss_mb");
    return metrics;
}

} // namespace image
