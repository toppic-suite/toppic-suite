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

#include "sql/ms_sql/ms_sql_writer.hpp"

#include <sqlite3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
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

// Compile a statement, aborting on failure (matching sql_util::execSql).
sqlite3_stmt* prepare(sqlite3* db, const std::string& sql) {
  sqlite3_stmt* stmt = nullptr;
  if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) != SQLITE_OK) {
    LOG_ERROR("Failed to prepare SQL: " << sql);
    LOG_ERROR("SQL error: " << sqlite3_errmsg(db));
    exit(EXIT_FAILURE);
  }
  return stmt;
}

// Run a fully-bound statement and reset it so the handle can be reused (every
// parameter is rebound before the next step, so no clear_bindings is needed).
void stepAndReset(sqlite3* db, sqlite3_stmt* stmt) {
  if (sqlite3_step(stmt) != SQLITE_DONE) {
    LOG_ERROR("SQL error: " << sqlite3_errmsg(db));
    exit(EXIT_FAILURE);
  }
  sqlite3_reset(stmt);
}

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

// Map the intensities to COLOR_NUM buckets on a log scale between the map's
// minimum and maximum intensity.
void setColors(std::vector<MsSqlPeak>& peaks, const MsSqlRange& range) {
  if (range.int_min <= 0.0 || range.int_max <= range.int_min) {
    for (MsSqlPeak& peak : peaks) {
      peak.color = 0;
    }
    return;
  }
  double log_min = std::log(range.int_min);
  double log_span = std::log(range.int_max) - log_min;
  for (MsSqlPeak& peak : peaks) {
    int idx = static_cast<int>(MsSqlWriter::COLOR_NUM *
                               (std::log(peak.inte) - log_min) / log_span);
    peak.color = std::clamp(idx, 0, MsSqlWriter::COLOR_NUM - 1);
  }
}

// Keep, for one retention-time row of the grid, the most intense peak of each
// m/z block. The kept peaks are appended to layer and extend layer_range.
void addRow(std::vector<MsSqlPeak>& row, double mz_min, double mz_size,
            std::vector<MsSqlPeak>& layer, MsSqlRange& layer_range) {
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
        extendRange(layer_range, peak);
        block_max_inte = peak.inte;
      }
    } else {
      layer.push_back(peak);
      extendRange(layer_range, peak);
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
// block. peaks is reordered by retention time.
std::vector<MsSqlPeak> downsample(std::vector<MsSqlPeak>& peaks,
                                  const MsSqlRange& range, double mz_size,
                                  double rt_size, MsSqlRange& layer_range) {
  std::sort(peaks.begin(), peaks.end(),
            [](const MsSqlPeak& a, const MsSqlPeak& b) { return a.rt < b.rt; });
  layer_range = emptyRange();
  layer_range.count = 0;
  std::vector<MsSqlPeak> layer;
  std::vector<MsSqlPeak> row;
  int row_num = 1;
  for (const MsSqlPeak& peak : peaks) {
    if (peak.rt >= range.rt_min + rt_size * row_num) {
      addRow(row, range.mz_min, mz_size, layer, layer_range);
      row.clear();
      row_num++;
    }
    row.push_back(peak);
  }
  addRow(row, range.mz_min, mz_size, layer, layer_range);
  // extendRange counted every kept peak, including the ones later replaced
  layer_range.count = static_cast<int>(layer.size());
  return layer;
}

// Order the peaks the way the (RETENTIONTIME, MZ) index is ordered, so that
// inserting them after the index exists only appends to it.
void sortByRtMz(std::vector<MsSqlPeak>& peaks) {
  std::sort(peaks.begin(), peaks.end(),
            [](const MsSqlPeak& a, const MsSqlPeak& b) {
              return a.rt < b.rt || (a.rt == b.rt && a.mz < b.mz);
            });
}

}  // namespace

void MsSqlWriter::write(std::vector<MsSqlPeak> peaks, int ms1_scan_num,
                        double mz_size, double rt_divider) {
  MsSqlRange range = computeRange(peaks);
  setColors(peaks, range);
  mz_size = std::min(mz_size, range.mz_max - range.mz_min);
  double rt_size = 0.0;
  if (ms1_scan_num > 0 && rt_divider > 0.0) {
    rt_size = static_cast<double>(range.rt_max - range.rt_min) /
              ms1_scan_num / rt_divider;
  }

  sql_util::execSql(db_, "BEGIN;");
  createConfigTable();
  int layer = 0;
  sortByRtMz(peaks);
  createLayerTable(layer);
  insertLayerPeaks(peaks, layer);
  insertConfig(range);
  std::cout << peaks.size() << " peaks written to PEAKS0." << std::endl;

  while (static_cast<int>(peaks.size()) >= MIN_LAYER_PEAKS && mz_size > 0.0 &&
         rt_size > 0.0) {
    MsSqlRange layer_range;
    std::vector<MsSqlPeak> layer_peaks =
        downsample(peaks, range, mz_size, rt_size, layer_range);
    mz_size *= MZ_SCALE;
    rt_size *= RT_SCALE;
    if (layer_peaks.size() == peaks.size()) {
      // no block held two peaks, so the layer would repeat the previous one
      continue;
    }
    layer++;
    sortByRtMz(layer_peaks);
    createLayerTable(layer);
    insertLayerPeaks(layer_peaks, layer);
    insertConfig(layer_range);
    std::cout << layer_peaks.size() << " peaks written to PEAKS" << layer << "."
              << std::endl;
    peaks = std::move(layer_peaks);
  }
  sql_util::execSql(db_, "COMMIT;");
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
  sqlite3_stmt* stmt =
      prepare(db_,
              "INSERT INTO CONFIG(MZMIN, MZMAX, RTMIN, RTMAX, INTMIN, INTMAX, "
              "COUNT) VALUES (?, ?, ?, ?, ?, ?, ?);");
  sqlite3_bind_double(stmt, 1, range.mz_min);
  sqlite3_bind_double(stmt, 2, range.mz_max);
  sqlite3_bind_int(stmt, 3, range.rt_min);
  sqlite3_bind_int(stmt, 4, range.rt_max);
  sqlite3_bind_double(stmt, 5, range.int_min);
  sqlite3_bind_double(stmt, 6, range.int_max);
  sqlite3_bind_int(stmt, 7, range.count);
  stepAndReset(db_, stmt);
  sqlite3_finalize(stmt);
}

void MsSqlWriter::createLayerTable(int layer) {
  std::string num = std::to_string(layer);
  sql_util::execSql(db_, "CREATE TABLE PEAKS" + num +
                             "(MZ REAL NOT NULL,"
                             "INTENSITY REAL NOT NULL,"
                             "RETENTIONTIME INT NOT NULL,"
                             "COLOR TINYINT NOT NULL);");
  // The index is created before the rows: they arrive in index order (see
  // sortByRtMz), so building it incrementally is cheaper than sorting the
  // whole table afterwards.
  sql_util::execSql(db_, "CREATE INDEX rtmz_index" + num + " ON PEAKS" + num +
                             " (RETENTIONTIME, MZ);");
}

void MsSqlWriter::insertLayerPeaks(const std::vector<MsSqlPeak>& peaks,
                                   int layer) {
  sqlite3_stmt* stmt =
      prepare(db_, "INSERT INTO PEAKS" + std::to_string(layer) +
                       "(MZ, INTENSITY, RETENTIONTIME, COLOR) "
                       "VALUES (?, ?, ?, ?);");
  for (const MsSqlPeak& peak : peaks) {
    sqlite3_bind_double(stmt, 1, peak.mz);
    sqlite3_bind_double(stmt, 2, peak.inte);
    sqlite3_bind_int(stmt, 3, peak.rt);
    sqlite3_bind_int(stmt, 4, peak.color);
    stepAndReset(db_, stmt);
  }
  sqlite3_finalize(stmt);
}

}  // namespace toppic
