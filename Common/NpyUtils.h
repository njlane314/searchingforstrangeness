#ifndef COMMON_NPYUTILS_H
#define COMMON_NPYUTILS_H

#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "art/Utilities/Exception.h"

namespace common {

inline void save_npy_f32_1d(const std::string &path,
                            const std::vector<float> &data) {
  const char magic[] = "\x93NUMPY";
  const uint8_t major = 1, minor = 0;

  std::ostringstream dict;
  dict << "{'descr': '<f4', 'fortran_order': False, 'shape': (" << data.size()
       << ",), }";
  std::string header = dict.str();

  size_t pre = sizeof(magic) - 1 + 2 + 2;
  size_t pad_len = 16 - ((pre + header.size()) % 16);
  if (pad_len == 16)
    pad_len = 0;
  header.append(pad_len, ' ');
  header.push_back('\n');

  uint16_t hlen = static_cast<uint16_t>(header.size());

  std::ofstream ofs(path, std::ios::binary);
  if (!ofs)
    throw art::Exception(art::errors::LogicError)
        << "Cannot open " << path << " for writing";

  ofs.write(magic, sizeof(magic) - 1);
  ofs.put(static_cast<char>(major));
  ofs.put(static_cast<char>(minor));
  ofs.write(reinterpret_cast<const char *>(&hlen), sizeof(hlen));
  ofs.write(header.data(), header.size());
  ofs.write(reinterpret_cast<const char *>(data.data()),
            static_cast<std::streamsize>(data.size() * sizeof(float)));
  if (!ofs)
    throw art::Exception(art::errors::LogicError) << "Short write to " << path;
}

inline void save_npy_f32_2d(const std::string &path,
                            const std::vector<std::vector<float>> &data) {
  if (data.empty())
    throw art::Exception(art::errors::LogicError)
        << "Cannot write empty array to " << path;

  size_t rows = data.size();
  size_t cols = data.front().size();
  for (const auto &row : data) {
    if (row.size() != cols)
      throw art::Exception(art::errors::LogicError)
          << "Inconsistent row size writing " << path;
  }

  const char magic[] = "\x93NUMPY";
  const uint8_t major = 1, minor = 0;

  std::ostringstream dict;
  dict << "{'descr': '<f4', 'fortran_order': False, 'shape': (" << rows << ", "
       << cols << ",), }";
  std::string header = dict.str();

  size_t pre = sizeof(magic) - 1 + 2 + 2;
  size_t pad_len = 16 - ((pre + header.size()) % 16);
  if (pad_len == 16)
    pad_len = 0;
  header.append(pad_len, ' ');
  header.push_back('\n');

  uint16_t hlen = static_cast<uint16_t>(header.size());

  std::ofstream ofs(path, std::ios::binary);
  if (!ofs)
    throw art::Exception(art::errors::LogicError)
        << "Cannot open " << path << " for writing";

  ofs.write(magic, sizeof(magic) - 1);
  ofs.put(static_cast<char>(major));
  ofs.put(static_cast<char>(minor));
  ofs.write(reinterpret_cast<const char *>(&hlen), sizeof(hlen));
  ofs.write(header.data(), header.size());

  for (const auto &row : data) {
    ofs.write(reinterpret_cast<const char *>(row.data()),
              static_cast<std::streamsize>(row.size() * sizeof(float)));
  }
  if (!ofs)
    throw art::Exception(art::errors::LogicError) << "Short write to " << path;
}

} // namespace common

#endif // COMMON_NPYUTILS_H
