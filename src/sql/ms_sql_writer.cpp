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

#include "sql/ms_sql_writer.hpp"

#include <sqlite3.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "common/util/logger.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace {

MsSqlRange emptyRange() {
  MsSqlRange range;
  range.mz_min = std::numeric_limits<double>::max();
  range.mz_max = std::numeric_limits<double>::lowest();
  range.rt_min = std::numeric_limits<int>::max();
  range.rt_max = std::numeric_limits<int>::lowest();
  range.int_min = std::numeric_limits<double>::max();
  range.int_max = std::numeric_limits<double>::lowest();
  return range;
}

void extendRange(MsSqlRange& range, const MsSqlPeak& peak) {
  range.mz_min = std::min(range.mz_min, peak.mz);
  range.mz_max = std::max(range.mz_max, peak.mz);
  range.rt_min = std::min(range.rt_min, peak.rt);
  range.rt_max = std::max(range.rt_max, peak.rt);
  range.int_min = std::min(range.int_min, peak.inte);
  range.int_max = std::max(range.int_max, peak.inte);
  range.count++;
}

MsSqlRange computeRange(const std::vector<MsSqlPeak>& peaks) {
  if (peaks.empty()) {
    return MsSqlRange();
  }
  MsSqlRange range = emptyRange();
  for (const MsSqlPeak& peak : peaks) {
    extendRange(range, peak);
  }
  return range;
}

// Keep, for one retention-time row of the grid, the most intense peak of each
// m/z block. The kept peaks are appended to layer.
void addRow(std::vector<MsSqlPeak>& row, double mz_min, double mz_size,
            std::vector<MsSqlPeak>& layer) {
  std::sort(row.begin(), row.end(),
            [](const MsSqlPeak& a, const MsSqlPeak& b) { return a.mz < b.mz; });
  double block_end = mz_min;  // upper m/z bound of the current block
  double block_max_inte = 0.0;
  bool first = true;
  for (const MsSqlPeak& peak : row) {
    if (peak.mz <= block_end) {
      // same block as the previous peak: keep only the more intense one
      if (peak.inte > block_max_inte) {
        if (!first && !layer.empty()) {
          layer.pop_back();
        }
        layer.push_back(peak);
        block_max_inte = peak.inte;
      }
    } else {
      layer.push_back(peak);
      block_max_inte = peak.inte;
      while (block_end < peak.mz) {
        block_end += mz_size;
      }
    }
    first = false;
  }
}

// One down-sampling step: cut the map into rt_size x mz_size blocks (from the
// map's minimum retention time and m/z) and keep the most intense peak of each
// block. peaks is reordered by retention time. layer_range is the range of the
// kept peaks, computed once the layer is complete: a peak replaced by a more
// intense one of its block must not widen it.
std::vector<MsSqlPeak> downsample(std::vector<MsSqlPeak>& peaks,
                                  const MsSqlRange& range, double mz_size,
                                  double rt_size, MsSqlRange& layer_range) {
  std::sort(peaks.begin(), peaks.end(),
            [](const MsSqlPeak& a, const MsSqlPeak& b) { return a.rt < b.rt; });
  std::vector<MsSqlPeak> layer;
  std::vector<MsSqlPeak> row;
  int row_num = 1;
  for (const MsSqlPeak& peak : peaks) {
    if (peak.rt >= range.rt_min + rt_size * row_num) {
      addRow(row, range.mz_min, mz_size, layer);
      row.clear();
      row_num++;
    }
    row.push_back(peak);
  }
  addRow(row, range.mz_min, mz_size, layer);
  layer_range = computeRange(layer);
  return layer;
}

// Order the peaks by retention time and m/z: the order the levels were
// written in as the tables PEAKS<k> with an index on (RETENTIONTIME, MZ)
// before toppic_test experiment 165. The down-sampling of the next level
// starts from it, and its sorts are not stable, so this order decides
// between peaks they find equal: kept, the levels are those of before.
void sortByRtMz(std::vector<MsSqlPeak>& peaks) {
  std::sort(peaks.begin(), peaks.end(),
            [](const MsSqlPeak& a, const MsSqlPeak& b) {
              return a.rt < b.rt || (a.rt == b.rt && a.mz < b.mz);
            });
}

// Put the peaks back in their input order (pos) in place: every swap moves
// one peak to its place.
void restoreInputOrder(std::vector<MsSqlPeak>& peaks) {
  for (std::size_t i = 0; i < peaks.size(); i++) {
    while (peaks[i].pos != static_cast<int>(i)) {
      std::size_t j = static_cast<std::size_t>(peaks[i].pos);
      std::swap(peaks[i], peaks[j]);
    }
  }
}

}  // namespace

void MsSqlWriter::write(std::vector<MsSqlPeak> peaks, int ms1_scan_num,
                        double mz_size, double rt_divider) {
  // The peaks of intensity <= 0 are not in the map (LEVEL -1). The others
  // keep their order, the order the levels have always been built from.
  std::vector<MsSqlPeak> unmapped;
  std::size_t peak_num = 0;
  for (const MsSqlPeak& peak : peaks) {
    if (peak.inte > 0.0) {
      peaks[peak_num++] = peak;
    } else {
      unmapped.push_back(peak);
    }
  }
  peaks.resize(peak_num);
  for (std::size_t i = 0; i < peak_num; i++) {
    peaks[i].pos = static_cast<int>(i);
  }
  // the highest level that keeps a peak, by its position in the input
  std::vector<std::int8_t> levels(peak_num, 0);

  MsSqlRange range = computeRange(peaks);
  mz_size = std::min(mz_size, range.mz_max - range.mz_min);
  double rt_size = 0.0;
  if (ms1_scan_num > 0 && rt_divider > 0.0) {
    rt_size = static_cast<double>(range.rt_max - range.rt_min) /
              ms1_scan_num / rt_divider;
  }

  sql_util::execSql(db_, "BEGIN;");
  createConfigTable();
  int level = 0;
  sortByRtMz(peaks);
  insertConfig(range);

  // each level is down-sampled from the one below (last); level 0 is kept
  // for ms1_peak
  std::vector<MsSqlPeak> last_level;
  std::vector<MsSqlPeak>* last = &peaks;
  while (static_cast<int>(last->size()) >= MIN_LAYER_PEAKS && mz_size > 0.0 &&
         rt_size > 0.0) {
    MsSqlRange layer_range;
    std::vector<MsSqlPeak> layer_peaks =
        downsample(*last, range, mz_size, rt_size, layer_range);
    mz_size *= MZ_SCALE;
    rt_size *= RT_SCALE;
    if (layer_peaks.size() == last->size()) {
      // no block held two peaks, so the level would repeat the previous one
      continue;
    }
    level++;
    sortByRtMz(layer_peaks);
    for (const MsSqlPeak& peak : layer_peaks) {
      levels[peak.pos] = static_cast<std::int8_t>(level);
    }
    insertConfig(layer_range);
    last_level = std::move(layer_peaks);
    last = &last_level;
  }
  std::vector<MsSqlPeak>().swap(last_level);

  // ms1_peak in key order (LEVEL, spec_id, peak_id): the input order within
  // a level
  restoreInputOrder(peaks);
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db_,
      "INSERT INTO ms1_peak(LEVEL, spec_id, peak_id, mz, intensity) "
      "VALUES (?, ?, ?, ?, ?);");
  auto insert = [&](int peak_level, const MsSqlPeak& peak) {
    sqlite3_bind_int(stmt, 1, peak_level);
    sqlite3_bind_int(stmt, 2, peak.spec_id);
    sqlite3_bind_int(stmt, 3, peak.peak_id);
    sqlite3_bind_double(stmt, 4, peak.mz);
    sqlite3_bind_double(stmt, 5, peak.inte);
    sql_util::stepAndReset(db_, stmt);
  };
  for (const MsSqlPeak& peak : unmapped) {
    insert(-1, peak);
  }
  for (int k = 0; k <= level; k++) {
    for (std::size_t i = 0; i < peak_num; i++) {
      if (levels[i] == k) {
        insert(k, peaks[i]);
      }
    }
  }
  sqlite3_finalize(stmt);
  sql_util::execSql(db_, "COMMIT;");
  std::cout << unmapped.size() + peak_num
            << " MS1 peaks written to ms1_peak with their 3D levels."
            << std::endl;
}

void MsSqlWriter::createConfigTable() {
  sql_util::execSql(db_,
                    "CREATE TABLE CONFIG("
                    "MZMIN REAL NOT NULL,"
                    "MZMAX REAL NOT NULL,"
                    "RTMIN INT NOT NULL,"
                    "RTMAX INT NOT NULL,"
                    "INTMIN REAL NOT NULL,"
                    "INTMAX REAL NOT NULL,"
                    "COUNT INT NOT NULL);");
}

void MsSqlWriter::insertConfig(const MsSqlRange& range) {
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db_,
      "INSERT INTO CONFIG(MZMIN, MZMAX, RTMIN, RTMAX, INTMIN, INTMAX, "
      "COUNT) VALUES (?, ?, ?, ?, ?, ?, ?);");
  sqlite3_bind_double(stmt, 1, range.mz_min);
  sqlite3_bind_double(stmt, 2, range.mz_max);
  sqlite3_bind_int(stmt, 3, range.rt_min);
  sqlite3_bind_int(stmt, 4, range.rt_max);
  sqlite3_bind_double(stmt, 5, range.int_min);
  sqlite3_bind_double(stmt, 6, range.int_max);
  sqlite3_bind_int(stmt, 7, range.count);
  sql_util::stepAndReset(db_, stmt);
  sqlite3_finalize(stmt);
}

}  // namespace toppic
