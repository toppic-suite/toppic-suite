// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "filter/massmatch/mass_match.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#ifndef _WIN32
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <zstd.h>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"

namespace toppic {


/* proteo_minus_water_masses[i]: the mass (minus water) of the ith proteoform
 * mass_2d[i]: a vector containing prefix residue masses of the ith proteoform
 * float_shift_2d[i]: a vector containing all possible shifts of the ith
 * proteoform pos_2d[i]: a vector containing the first residue position for each
 * shift.
 */
MassMatch::MassMatch(std::vector<double>& proteo_minus_water_masses,
                     std::vector<std::vector<int>>& mass_2d,
                     std::vector<std::vector<double>>& float_shift_2d,
                     std::vector<std::vector<int>>& pos_2d,
                     double max_proteoform_mass, double scale) {
  scale_ = scale;
  LOG_DEBUG("Scale: " << scale_);

  proteo_num_ = proteo_minus_water_masses.size();
  LOG_DEBUG("Proteoform number: " << proteo_minus_water_masses.size());
  proteo_minus_water_masses_ = proteo_minus_water_masses;

  col_num_ = max_proteoform_mass * scale_;
  LOG_DEBUG("column number: " << col_num_);

  LOG_DEBUG("start init");
  initProteoformBeginEnds(float_shift_2d);
  LOG_DEBUG("row number: " << row_num_);

  LOG_DEBUG("init indexes");
  initIndexes(mass_2d, float_shift_2d, pos_2d);
  initColIndexBounds();
}

int MassMatch::write_version_ = 2;

namespace {

// The version 2 index file (toppic_test experiment 150): the magic string,
// the scalars, then the seven arrays, each as its element count, its
// zstd frame size and the frame. The column bounds are stored once (col_num
// + 1 entries) instead of the begins and ends of version 1 (80 MB of the
// 400 MB of a zero-shift diagonal file of the human database); the row ids
// compress by about 1.15 and the bounds by about 30.
const char kIndexMagicV2[8] = {'T', 'P', 'M', 'M', 'I', 'D', 'X', '2'};
const int kZstdLevel = 1;

void writeRaw(std::ofstream& out, const void* data, size_t size) {
  out.write(static_cast<const char*>(data), static_cast<std::streamsize>(size));
}

template <typename T>
void writeArrayV2(std::ofstream& out, const std::vector<T>& vec,
                  std::vector<char>& buffer) {
  uint64_t count = vec.size();
  size_t raw_size = count * sizeof(T);
  size_t bound = ZSTD_compressBound(raw_size);
  if (buffer.size() < bound) {
    buffer.resize(bound);
  }
  size_t comp_size = 0;
  if (raw_size > 0) {
    comp_size = ZSTD_compress(buffer.data(), bound, vec.data(), raw_size,
                              kZstdLevel);
    if (ZSTD_isError(comp_size)) {
      LOG_ERROR("Index compression failed: " << ZSTD_getErrorName(comp_size));
      exit(EXIT_FAILURE);
    }
  }
  uint64_t comp_size64 = comp_size;
  writeRaw(out, &count, sizeof(count));
  writeRaw(out, &comp_size64, sizeof(comp_size64));
  if (comp_size > 0) {
    writeRaw(out, buffer.data(), comp_size);
  }
}

template <typename T>
void readArrayV2(const char* data, size_t data_size, size_t& pos,
                 std::vector<T>& vec, const std::string& file_path) {
  if (pos + 16 > data_size) {
    LOG_ERROR("Truncated index file: " << file_path);
    exit(EXIT_FAILURE);
  }
  uint64_t count, comp_size;
  std::memcpy(&count, data + pos, 8);
  std::memcpy(&comp_size, data + pos + 8, 8);
  pos += 16;
  if (pos + comp_size > data_size) {
    LOG_ERROR("Truncated index file: " << file_path);
    exit(EXIT_FAILURE);
  }
  vec.resize(count);
  if (count > 0) {
    size_t raw_size = count * sizeof(T);
    size_t got = ZSTD_decompress(vec.data(), raw_size, data + pos, comp_size);
    if (ZSTD_isError(got) || got != raw_size) {
      LOG_ERROR("Index decompression failed: " << file_path);
      exit(EXIT_FAILURE);
    }
  }
  pos += comp_size;
}

// Reads the whole file into buffer with read(2) (not the C++ stream, so that
// the page cache can be told to drop the file afterwards) and returns the
// bytes read; the index files stream through the page cache once per stage,
// 100 GB and more for the human database, and without the release the
// kernel evicts under pressure for the whole stage (free memory at 0.3 GB,
// 13-19% system time in toppic_test experiment 99).
size_t readWholeFile(const std::string& file_path, std::vector<char>& buffer) {
#ifndef _WIN32
  int fd = ::open(file_path.c_str(), O_RDONLY);
  if (fd < 0) {
    return 0;
  }
  struct stat st;
  if (::fstat(fd, &st) != 0) {
    ::close(fd);
    return 0;
  }
  size_t size = static_cast<size_t>(st.st_size);
  if (buffer.size() < size) {
    buffer.resize(size);
  }
  size_t pos = 0;
  while (pos < size) {
    ssize_t got = ::read(fd, buffer.data() + pos, size - pos);
    if (got <= 0) {
      break;
    }
    pos += static_cast<size_t>(got);
  }
#ifdef POSIX_FADV_DONTNEED
  ::posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
#endif
  ::close(fd);
  return pos;
#else
  std::ifstream in(file_path, std::ifstream::binary | std::ifstream::ate);
  if (!in.is_open()) {
    return 0;
  }
  size_t size = static_cast<size_t>(in.tellg());
  in.seekg(0);
  if (buffer.size() < size) {
    buffer.resize(size);
  }
  in.read(buffer.data(), static_cast<std::streamsize>(size));
  return static_cast<size_t>(in.gcount());
#endif
}

}  // namespace

void MassMatch::serializeMassMatch(std::string file_name,
                                   std::string dir_name) {
  std::string file_path = dir_name + file_util::getFileSeparator() + file_name;
  if (write_version_ == 2) {
    serializeMassMatchV2(file_path);
  } else {
    serializeMassMatchV1(file_path);
  }
}

void MassMatch::serializeMassMatchV1(const std::string& file_path) {
  std::ofstream new_file(file_path, std::ofstream::binary);

  if (new_file.is_open()) {
    // the begins and ends of the archive, from the bounds (initColIndexBounds
    // freed them)
    if (col_index_begins_.empty() && col_num_ > 0) {
      col_index_begins_.resize(col_num_);
      col_index_ends_.resize(col_num_);
      for (int i = 0; i < col_num_; i++) {
        col_index_begins_[i] = col_index_bounds_[i];
        col_index_ends_[i] = col_index_bounds_[i + 1] - 1;
      }
    }
    boost::archive::binary_oarchive oa(new_file, std::ios::binary);
    oa << scale_;
    oa << proteo_num_;
    oa << col_num_;
    oa << row_num_;

    oa << proteo_row_begins_;
    oa << proteo_row_ends_;
    oa << proteo_minus_water_masses_;
    oa << row_proteo_ids_;
    oa << trunc_shifts_;

    oa << col_index_begins_;
    oa << col_index_ends_;
    oa << col_indexes_;

    new_file.close();
    std::vector<int>().swap(col_index_begins_);
    std::vector<int>().swap(col_index_ends_);
  } else {
    LOG_ERROR("Failed to open the index file:" << file_path);
    exit(EXIT_FAILURE);
  }
}

void MassMatch::serializeMassMatchV2(const std::string& file_path) {
  std::ofstream out(file_path, std::ofstream::binary);
  if (!out.is_open()) {
    LOG_ERROR("Failed to open the index file:" << file_path);
    exit(EXIT_FAILURE);
  }
  writeRaw(out, kIndexMagicV2, sizeof(kIndexMagicV2));
  writeRaw(out, &scale_, sizeof(scale_));
  int32_t scalars[3] = {proteo_num_, col_num_, row_num_};
  writeRaw(out, scalars, sizeof(scalars));
  std::vector<char> buffer;
  writeArrayV2(out, proteo_row_begins_, buffer);
  writeArrayV2(out, proteo_row_ends_, buffer);
  writeArrayV2(out, proteo_minus_water_masses_, buffer);
  writeArrayV2(out, row_proteo_ids_, buffer);
  writeArrayV2(out, trunc_shifts_, buffer);
  writeArrayV2(out, col_index_bounds_, buffer);
  writeArrayV2(out, col_indexes_, buffer);
  out.close();
  if (!out) {
    LOG_ERROR("Failed to write the index file:" << file_path);
    exit(EXIT_FAILURE);
  }
}

void MassMatch::deserializeMassMatch(std::string new_file,
                                     std::string dir_name) {
  std::string file_path = dir_name + file_util::getFileSeparator() + new_file;
  // the format: the version 2 magic string, else the version 1 archive
  char head[sizeof(kIndexMagicV2)] = {0};
  {
    std::ifstream probe(file_path, std::ifstream::binary);
    if (!probe.is_open()) {
      LOG_ERROR("Failed to open the index file:" << file_path);
      exit(EXIT_FAILURE);
    }
    probe.read(head, sizeof(head));
  }
  if (std::memcmp(head, kIndexMagicV2, sizeof(kIndexMagicV2)) == 0) {
    deserializeMassMatchV2(file_path);
  } else {
    deserializeMassMatchV1(file_path);
  }
}

void MassMatch::deserializeMassMatchV1(const std::string& file_path) {
  std::ifstream file_to_read(file_path, std::ifstream::binary);

  if (file_to_read.is_open()) {
    boost::archive::binary_iarchive ia(file_to_read, std::ios::binary);

    ia >> scale_;
    ia >> proteo_num_;
    ia >> col_num_;
    ia >> row_num_;

    ia >> proteo_row_begins_;
    ia >> proteo_row_ends_;
    ia >> proteo_minus_water_masses_;
    ia >> row_proteo_ids_;
    ia >> trunc_shifts_;

    ia >> col_index_begins_;
    ia >> col_index_ends_;
    ia >> col_indexes_;

    file_to_read.close();
    initColIndexBounds();
  } else {
    LOG_ERROR("Failed to open the index file:" << file_path);
    exit(EXIT_FAILURE);
  }
}

void MassMatch::deserializeMassMatchV2(const std::string& file_path) {
  // one buffer per thread: the files of a stage are read one after another
  // by the same thread
  static thread_local std::vector<char> buffer;
  size_t size = readWholeFile(file_path, buffer);
  size_t header_size = sizeof(kIndexMagicV2) + sizeof(double) + 3 * 4;
  if (size < header_size) {
    LOG_ERROR("Failed to read the index file:" << file_path);
    exit(EXIT_FAILURE);
  }
  const char* data = buffer.data();
  size_t pos = sizeof(kIndexMagicV2);
  std::memcpy(&scale_, data + pos, sizeof(double));
  pos += sizeof(double);
  int32_t scalars[3];
  std::memcpy(scalars, data + pos, sizeof(scalars));
  pos += sizeof(scalars);
  proteo_num_ = scalars[0];
  col_num_ = scalars[1];
  row_num_ = scalars[2];
  readArrayV2(data, size, pos, proteo_row_begins_, file_path);
  readArrayV2(data, size, pos, proteo_row_ends_, file_path);
  readArrayV2(data, size, pos, proteo_minus_water_masses_, file_path);
  readArrayV2(data, size, pos, row_proteo_ids_, file_path);
  readArrayV2(data, size, pos, trunc_shifts_, file_path);
  readArrayV2(data, size, pos, col_index_bounds_, file_path);
  readArrayV2(data, size, pos, col_indexes_, file_path);
  if (static_cast<int>(col_index_bounds_.size()) != col_num_ + 1) {
    LOG_ERROR("Bad column bounds in the index file:" << file_path);
    exit(EXIT_FAILURE);
  }
}

void MassMatch::initProteoformBeginEnds(
    std::vector<std::vector<double>>& shift_2d) {
  // no need to init
  proteo_row_begins_.resize(proteo_num_);
  proteo_row_ends_.resize(proteo_num_);
  int pnt = 0;
  for (int i = 0; i < proteo_num_; i++) {
    proteo_row_begins_[i] = pnt;
    int len = shift_2d[i].size();
    proteo_row_ends_[i] = pnt + len - 1;
    pnt += len;
  }
  row_num_ = pnt;
  row_proteo_ids_.resize(row_num_);
  trunc_shifts_.resize(row_num_);
  for (int i = 0; i < proteo_num_; i++) {
    for (int j = proteo_row_begins_[i]; j <= proteo_row_ends_[i]; j++) {
      row_proteo_ids_[j] = i;
      int pos = j - proteo_row_begins_[i];
      trunc_shifts_[j] = shift_2d[i][pos];
    }
  }
}

inline std::vector<std::vector<int>> convertToInt(
    std::vector<std::vector<double>>& mass_2d, double scale) {
  std::vector<std::vector<int>> result;
  for (size_t i = 0; i < mass_2d.size(); i++) {
    std::vector<int> int_masses;
    for (size_t j = 0; j < mass_2d[i].size(); j++) {
      int value = std::floor(mass_2d[i][j] * scale + 0.5);
      int_masses.push_back(value);
    }
    result.push_back(int_masses);
  }
  return result;
}

void MassMatch::compColumnMatchNums(std::vector<std::vector<int>>& mass_2d,
                                    std::vector<std::vector<int>>& shift_2d,
                                    std::vector<std::vector<int>>& pos_2d,
                                    std::vector<int>& col_match_nums) {
  size_t proteo_num = mass_2d.size();
  for (size_t i = 0; i < proteo_num; i++) {
    for (size_t s = 0; s < shift_2d[i].size(); s++) {
      for (size_t cur = pos_2d[i][s]; cur < mass_2d[i].size(); cur++) {
        int shift_mass = mass_2d[i][cur] + shift_2d[i][s];
        if (shift_mass > 0) {
          if (shift_mass < col_num_) {
            col_match_nums[shift_mass]++;
          } else {
            break;
          }
        }
      }
    }
  }
}

void MassMatch::fillColumnIndex(std::vector<std::vector<int>>& mass_2d,
                                std::vector<std::vector<int>>& shift_2d,
                                std::vector<std::vector<int>>& pos_2d,
                                std::vector<int>& col_index_pnts) {
  for (size_t i = 0; i < mass_2d.size(); i++) {
    for (size_t s = 0; s < shift_2d[i].size(); s++) {
      for (size_t cur = pos_2d[i][s]; cur < mass_2d[i].size(); cur++) {
        int shift_mass = mass_2d[i][cur] + shift_2d[i][s];
        if (shift_mass > 0) {
          if (shift_mass < col_num_) {
            col_indexes_[col_index_pnts[shift_mass]] =
                proteo_row_begins_[i] + s;
            col_index_pnts[shift_mass]++;
          } else {
            break;
          }
        }
      }
    }
  }
}

// The column index built range by range. Row r = (proteoform i, shift s)
// contributes the entries (mass + shift, r) for its masses from the start
// position on while 0 < mass + shift < col_num, in the order of the masses;
// with nondecreasing masses those entries are one run of the row's masses,
// and the part of the run that falls into a range of columns [c0, c1)
// starts where the part for the previous range ended. The columns are
// therefore processed in ranges of kColRangeWidth: a count pass and a fill
// pass over the rows' runs in that range, whose counts and write pointers
// stay in cache, and the running total of the entries gives the begins
// and ends. The result is the same as two passes over all columns at once
// (compColumnMatchNums, fillColumnIndex), whose random accesses to arrays
// of col_num entries miss the cache on every entry.
namespace {
const int kColRangeWidth = 1 << 16;

void buildColumnIndexByRange(const std::vector<std::vector<int>>& mass_2d,
                             const std::vector<std::vector<int>>& shift_2d,
                             const std::vector<std::vector<int>>& pos_2d,
                             int col_num, int row_num,
                             const std::vector<int>& proteo_row_begins,
                             std::vector<int>& col_index_begins,
                             std::vector<int>& col_index_ends,
                             std::vector<int>& col_indexes) {
  // each row's run: its current position and its end
  std::vector<int> row_cur(row_num);
  std::vector<int> row_end(row_num);
  size_t entry_num = 0;
  for (size_t i = 0; i < mass_2d.size(); i++) {
    const std::vector<int>& masses = mass_2d[i];
    for (size_t s = 0; s < shift_2d[i].size(); s++) {
      int r = proteo_row_begins[i] + s;
      int shift = shift_2d[i][s];
      auto begin = masses.begin() + pos_2d[i][s];
      // mass + shift > 0 and mass + shift < col_num
      auto lo = std::lower_bound(begin, masses.end(), 1 - shift);
      auto hi = std::lower_bound(lo, masses.end(), col_num - shift);
      row_cur[r] = lo - masses.begin();
      row_end[r] = hi - masses.begin();
      entry_num += hi - lo;
    }
  }
  col_index_begins.resize(col_num);
  col_index_ends.resize(col_num);
  col_indexes.resize(entry_num);

  std::vector<int> counts(kColRangeWidth);
  std::vector<int> pnts(kColRangeWidth);
  std::vector<int> range_end(row_num);
  int base = 0;
  for (int c0 = 0; c0 < col_num; c0 += kColRangeWidth) {
    int c1 = std::min(c0 + kColRangeWidth, col_num);
    std::fill(counts.begin(), counts.begin() + (c1 - c0), 0);
    for (size_t i = 0; i < mass_2d.size(); i++) {
      const int* masses = mass_2d[i].data();
      for (size_t s = 0; s < shift_2d[i].size(); s++) {
        int r = proteo_row_begins[i] + s;
        int shift = shift_2d[i][s];
        int cur = row_cur[r];
        int end = row_end[r];
        while (cur < end && masses[cur] + shift < c1) {
          counts[masses[cur] + shift - c0]++;
          cur++;
        }
        range_end[r] = cur;
      }
    }
    for (int c = c0; c < c1; c++) {
      col_index_begins[c] = base;
      pnts[c - c0] = base;
      col_index_ends[c] = base + counts[c - c0] - 1;
      base += counts[c - c0];
    }
    for (size_t i = 0; i < mass_2d.size(); i++) {
      const int* masses = mass_2d[i].data();
      for (size_t s = 0; s < shift_2d[i].size(); s++) {
        int r = proteo_row_begins[i] + s;
        int shift = shift_2d[i][s];
        int end = range_end[r];
        for (int cur = row_cur[r]; cur < end; cur++) {
          col_indexes[pnts[masses[cur] + shift - c0]++] = r;
        }
        row_cur[r] = end;
      }
    }
  }
}
}  // namespace

void MassMatch::initIndexes(std::vector<std::vector<int>>& mass_2d,
                            std::vector<std::vector<double>>& float_shift_2d,
                            std::vector<std::vector<int>>& pos_2d) {
  std::vector<std::vector<int>> shift_2d = convertToInt(float_shift_2d, scale_);
  LOG_DEBUG("column num " << col_num_);
  rows_nondecreasing_ = true;
  for (size_t i = 0; i < mass_2d.size(); i++) {
    if (!std::is_sorted(mass_2d[i].begin(), mass_2d[i].end())) {
      rows_nondecreasing_ = false;
      break;
    }
  }
  if (rows_nondecreasing_) {
    buildColumnIndexByRange(mass_2d, shift_2d, pos_2d, col_num_, row_num_,
                            proteo_row_begins_, col_index_begins_,
                            col_index_ends_, col_indexes_);
    LOG_DEBUG("indexes size: " << col_indexes_.size());
    return;
  }
  // masses that are not nondecreasing (no index of the suite makes them):
  // the two passes over all columns
  std::vector<int> col_match_nums(col_num_, 0);
  // no need to initialize
  std::vector<int> col_index_pnts(col_num_);
  col_index_begins_.resize(col_num_);
  col_index_ends_.resize(col_num_);

  compColumnMatchNums(mass_2d, shift_2d, pos_2d, col_match_nums);

  int pnt = 0;
  for (int i = 0; i < col_num_; i++) {
    col_index_begins_[i] = pnt;
    col_index_pnts[i] = pnt;
    col_index_ends_[i] = pnt + col_match_nums[i] - 1;
    pnt += col_match_nums[i];
  }
  // no need to initialize
  col_indexes_.resize(pnt, 0);
  LOG_DEBUG("indexes size: " << pnt);
  fillColumnIndex(mass_2d, shift_2d, pos_2d, col_index_pnts);
}

void MassMatch::initColIndexBounds() {
  col_index_bounds_ = std::move(col_index_begins_);
  col_index_bounds_.push_back(static_cast<int>(col_indexes_.size()));
  // the ends are the next column's begins; the scoring reads the bounds
  std::vector<int>().swap(col_index_begins_);
  std::vector<int>().swap(col_index_ends_);
}

bool MassMatch::canTruncateColumns(double max_proteoform_mass) const {
  int col_num = max_proteoform_mass * scale_;
  return rows_nondecreasing_ && col_num >= 0 && col_num <= col_num_;
}

void MassMatch::truncateColumns(double max_proteoform_mass) {
  if (!canTruncateColumns(max_proteoform_mass)) {
    LOG_ERROR("The index cannot be truncated to the proteoform mass "
              << max_proteoform_mass);
    exit(EXIT_FAILURE);
  }
  // the same conversion as the constructor's
  int col_num = max_proteoform_mass * scale_;
  // the entries of the columns below col_num come first
  int entry_num = col_index_bounds_[col_num];
  col_index_bounds_.resize(col_num + 1);
  col_indexes_.resize(entry_num);
  col_num_ = col_num;
}

void MassMatch::compScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors,
    std::vector<short>& scores) {
  compScores(pref_mass_errors, 0, 0.0, scores);
}

void MassMatch::compScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors, int start,
    double shift, std::vector<short>& scores, int col_limit) {
  // The column bounds (40 MB per index) and the row ids (hundreds of MB)
  // are touched at the masses' positions only, so every window misses the
  // cache; the bounds of the window a few masses ahead and the first row
  // ids of the next window are prefetched, as the diagonal filter's
  // compDiagMaxScoresImpl does (toppic_test experiments 16 and 150).
  const int* bounds = col_index_bounds_.data();
  const int* entries = col_indexes_.data();
  short* s = scores.data();
  const int n = static_cast<int>(pref_mass_errors.size());
  const int col_num =
      col_limit >= 0 && col_limit < col_num_ ? col_limit : col_num_;
  const int pf_bound = 8;  // masses ahead for the column bound prefetch
  const int pf_entry = 3;  // masses ahead for the row id prefetch
  for (int i = start; i < n; i++) {
    if (i + pf_bound < n) {
      int m2 = pref_mass_errors[i + pf_bound].first + shift;
      int left2 = m2 - pref_mass_errors[i + pf_bound].second;
      int right2 = m2 + pref_mass_errors[i + pf_bound].second;
      if (left2 < 0) {
        left2 = 0;
      }
      if (right2 >= 0 && right2 < col_num) {
        __builtin_prefetch(bounds + left2);
        __builtin_prefetch(bounds + right2 + 1);
      }
    }
    if (i + pf_entry < n) {
      int m2 = pref_mass_errors[i + pf_entry].first + shift;
      int left2 = m2 - pref_mass_errors[i + pf_entry].second;
      int right2 = m2 + pref_mass_errors[i + pf_entry].second;
      if (left2 < 0) {
        left2 = 0;
      }
      if (right2 >= 0 && right2 < col_num) {
        int begin2 = bounds[left2];
        __builtin_prefetch(entries + begin2);
        __builtin_prefetch(entries + begin2 + 16);
      }
    }
    int m = pref_mass_errors[i].first + shift;
    // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
    int left = m - pref_mass_errors[i].second;
    if (left < 0) {
      left = 0;
    }

    int right = m + pref_mass_errors[i].second;

    if (right < 0 || right >= col_num) {
      continue;
    }

    int begin_index = bounds[left];
    int end_index = bounds[right + 1];

    for (int j = begin_index; j < end_index; j++) {
      s[entries[j]]++;
    }
  }
}

// The anchor loop of DiagFilter::compute for the anchors [anchor_begin,
// anchor_end) with the score array of the current anchor in Score. The rows
// whose score reaches the threshold are noted by the increment loop itself
// (a rarely taken branch) and only they are folded into max_scores after the
// anchor, and the array is cleared with memset, so the per-anchor work on
// the whole array is one sequential write instead of a fill, a read and a
// read-modify-write of 16-bit scores. The column bounds and the first
// entries of the next pairs are prefetched: each pair reads them from memory
// the hardware prefetcher cannot predict.
template <typename Score>
void MassMatch::compDiagMaxScoresImpl(
    const std::vector<std::pair<int, int>>& pref_mass_errors, int anchor_begin,
    int anchor_end, int threshold, std::vector<Score>& scores,
    std::vector<int>& hot_rows, std::vector<short>& max_scores) {
  scores.assign(row_num_, 0);
  const int* bounds = col_index_bounds_.data();
  const int* entries = col_indexes_.data();
  Score* s = scores.data();
  short* mx = max_scores.data();
  const int n = static_cast<int>(pref_mass_errors.size());
  const int col_num = col_num_;
  const int row_num = row_num_;
  const int pf_bound = 8;  // pairs ahead for the column bound prefetch
  const int pf_entry = 3;  // pairs ahead for the entry prefetch
  const Score t = static_cast<Score>(threshold);
  for (int i = anchor_begin; i < anchor_end; i++) {
    if (i > anchor_begin) {
      std::memset(s, 0, sizeof(Score) * row_num);
    }
    hot_rows.clear();
    int shift = -pref_mass_errors[i].first;
    for (int k = i; k < n; k++) {
      if (k + pf_bound < n) {
        int m2 = pref_mass_errors[k + pf_bound].first + shift;
        int left2 = m2 - pref_mass_errors[k + pf_bound].second;
        int right2 = m2 + pref_mass_errors[k + pf_bound].second;
        if (left2 < 0) {
          left2 = 0;
        }
        if (right2 >= 0 && right2 < col_num) {
          __builtin_prefetch(bounds + left2);
          __builtin_prefetch(bounds + right2 + 1);
        }
      }
      if (k + pf_entry < n) {
        int m2 = pref_mass_errors[k + pf_entry].first + shift;
        int left2 = m2 - pref_mass_errors[k + pf_entry].second;
        int right2 = m2 + pref_mass_errors[k + pf_entry].second;
        if (left2 < 0) {
          left2 = 0;
        }
        if (right2 >= 0 && right2 < col_num) {
          int begin2 = bounds[left2];
          __builtin_prefetch(entries + begin2);
          __builtin_prefetch(entries + begin2 + 16);
        }
      }
      int m = pref_mass_errors[k].first + shift;
      int left = m - pref_mass_errors[k].second;
      if (left < 0) {
        left = 0;
      }
      int right = m + pref_mass_errors[k].second;
      if (right < 0 || right >= col_num) {
        continue;
      }
      int begin_index = bounds[left];
      int end_index = bounds[right + 1];
      for (int j = begin_index; j < end_index; j++) {
        int row = entries[j];
        if (++s[row] == t) {
          hot_rows.push_back(row);
        }
      }
    }
    for (size_t h = 0; h < hot_rows.size(); h++) {
      int row = hot_rows[h];
      if (s[row] > mx[row]) {
        mx[row] = static_cast<short>(s[row]);
      }
    }
  }
}

// An anchor's score for a row is at most the number of pairs of the anchor
// that hit the row, and a row is hit at most once by a pair when the column
// windows are narrower than the smallest residue mass (a row's entries are
// at least one residue mass apart). So the anchors with at most 255 later
// peaks are scored with 8-bit arrays and the earlier ones with 16-bit arrays.
void MassMatch::compDiagMaxScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors, int threshold,
    std::vector<short>& max_scores) {
  // scratch arrays, kept per thread across spectrum sets and blocks
  static thread_local std::vector<uint8_t> scores8;
  static thread_local std::vector<uint16_t> scores16;
  static thread_local std::vector<int> hot_rows;
  max_scores.assign(row_num_, 0);
  const int n = static_cast<int>(pref_mass_errors.size());
  int max_error = 0;
  for (int k = 0; k < n; k++) {
    max_error = std::max(max_error, pref_mass_errors[k].second);
  }
  // glycine, 57.02 Da, scaled by the index scale
  const int min_residue_mass = static_cast<int>(57.0 * scale_);
  int first8 = n;
  if (2 * max_error + 1 < min_residue_mass) {
    first8 = std::max(0, n - 255);
  }
  if (first8 > 0) {
    compDiagMaxScoresImpl(pref_mass_errors, 0, first8, threshold, scores16,
                          hot_rows, max_scores);
  }
  if (first8 < n) {
    compDiagMaxScoresImpl(pref_mass_errors, first8, n, threshold, scores8,
                          hot_rows, max_scores);
  }
}

template <bool TrackPrecRows>
void MassMatch::updatePrecScoreImpl(const std::pair<int, int>& mass_error,
                                    std::vector<short>& scores,
                                    std::vector<int>* prec_rows) {
  // precursor mass
  int begin_index, end_index;
  int m = mass_error.first;
  // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
  int left = m - mass_error.second;
  if (left < 0) {
    left = 0;
  }

  int right = m + mass_error.second;
  if (right >= col_num_) {
    right = col_num_ - 1;
  }

  if (left <= right) {
    // update scores
    begin_index = col_index_bounds_[left];
    end_index = col_index_bounds_[right + 1] - 1;
    for (int j = begin_index; j <= end_index; j++) {
      int row = col_indexes_[j];
      if (scores[row] < getPrecursorMatchScore()) {
        scores[row] += getPrecursorMatchScore();
        if (TrackPrecRows) {
          prec_rows->push_back(row);
        }
      }
    }
  }
}

void MassMatch::updatePrecScore(const std::pair<int, int> mass_error,
                                std::vector<short>& scores) {
  updatePrecScoreImpl<false>(mass_error, scores, nullptr);
}

void MassMatch::compMatchScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors,
    const std::pair<int, int>& prec_minus_water_mass_error,
    std::vector<short>& scores) {
  compScores(pref_mass_errors, 0, 0.0, scores);
  updatePrecScore(prec_minus_water_mass_error, scores);
}

void MassMatch::compMatchScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors,
    const std::vector<std::pair<int, int>>& prec_minus_water_mass_errors,
    std::vector<short>& scores) {
  compScores(pref_mass_errors, 0, 0.0, scores);
  // precursor mass
  for (size_t i = 0; i < prec_minus_water_mass_errors.size(); i++) {
    std::pair<int, int> mass_error = prec_minus_water_mass_errors[i];
    updatePrecScore(mass_error, scores);
  }
}

void MassMatch::compMatchScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors,
    const std::vector<std::pair<int, int>>& prec_minus_water_mass_errors,
    std::vector<short>& scores, std::vector<int>& prec_rows) {
  compScores(pref_mass_errors, 0, 0.0, scores);
  compPrecScores(prec_minus_water_mass_errors, scores, prec_rows);
}

void MassMatch::compPrecScores(
    const std::vector<std::pair<int, int>>& prec_minus_water_mass_errors,
    std::vector<short>& scores, std::vector<int>& prec_rows) {
  for (size_t i = 0; i < prec_minus_water_mass_errors.size(); i++) {
    updatePrecScoreImpl<true>(prec_minus_water_mass_errors[i], scores,
                              &prec_rows);
  }
}

} /* namespace toppic */
