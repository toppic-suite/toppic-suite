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

#ifndef TOPPIC_SQL_MS_SQL_WRITER_HPP_
#define TOPPIC_SQL_MS_SQL_WRITER_HPP_

#include <sqlite3.h>

#include <string>
#include <vector>

namespace toppic {

// One MS1 peak of an LC-MS map.
struct MsSqlPeak {
  double mz = 0.0;
  double inte = 0.0;
  int rt = 0;       // retention time in milliseconds
  int spec_id = 0;  // the peak's key in ms1_peak: its spectrum's id and
  int peak_id = 0;  // its index in the spectrum's peak list
  int pos = 0;      // set by MsSqlWriter::write: the position in its input
};

// The m/z, retention-time and intensity ranges of a set of peaks and their
// count: one row of the CONFIG table.
struct MsSqlRange {
  double mz_min = 0.0;
  double mz_max = 0.0;
  int rt_min = 0;
  int rt_max = 0;
  double int_min = 0.0;
  double int_max = 0.0;
  int count = 0;
};

// Writes the MS1 peaks of an LC-MS map into the SQLite database with the
// levels read by the 3D visualization (topfd --sql-3d; toppic_test
// experiment 165): the table ms1_peak and its CONFIG table.
//
// Level 0 of the map is every MS1 peak of intensity > 0, and levels 1, 2,
// ... are progressively down-sampled copies of it. The map is cut into a
// grid of rt_size x mz_size blocks and only the most intense peak of each
// block is kept; both block sizes double for each further level, until a
// level has fewer than MIN_LAYER_PEAKS peaks. The levels are nested, so
// ms1_peak stores each peak once with LEVEL, the highest level that keeps
// it (-1 for a peak of intensity <= 0, outside the map): level k is the rows
// with LEVEL >= k. Row k + 1 of the CONFIG table is the range of level k.
// The table is keyed on (LEVEL, spec_id, peak_id) (WITHOUT ROWID, created
// by TopfdPara::createSqlDb), so the peaks of level k in a range of spectra
// are one key range per LEVEL >= k, and a spectrum's peaks one key range per
// LEVEL. A peak's retention time is its spectrum's; the colour of the 3D
// view is computed by the reader from the intensity and CONFIG's first row.
// (Before experiment 165, level k was a table PEAKS<k> with the retention
// time and a COLOR column, and ms1_peak a separate table without LEVEL.)
class MsSqlWriter {
 public:
  // The sqlite3 connection is owned by the caller, which opens and closes it;
  // it must not hold an open transaction while write() runs.
  explicit MsSqlWriter(sqlite3* db) : db_(db) {}

  MsSqlWriter(const MsSqlWriter&) = delete;
  MsSqlWriter& operator=(const MsSqlWriter&) = delete;

  // Writes CONFIG and ms1_peak (which must be empty) in one transaction.
  // peaks is every MS1 peak, in the order of the ms1_peak key (spec_id,
  // peak_id) so that the rows of a level are appended. mz_size is the m/z
  // width of a grid block of the first down-sampled level (capped at the
  // map's m/z range); the block's retention-time height is the mean MS1 scan
  // interval divided by rt_divider.
  void write(std::vector<MsSqlPeak> peaks, int ms1_scan_num, double mz_size,
             double rt_divider);

  static constexpr double DEFAULT_MZ_SIZE = 0.05;
  static constexpr double DEFAULT_RT_DIVIDER = 1.0;

 private:
  void createConfigTable();
  void insertConfig(const MsSqlRange& range);

  sqlite3* db_;  // not owned

  // A level with fewer peaks than this is the last one.
  static constexpr int MIN_LAYER_PEAKS = 3000;
  // Growth of the grid blocks from one level to the next.
  static constexpr double MZ_SCALE = 2.0;
  static constexpr double RT_SCALE = 2.0;
};

}  // namespace toppic

#endif
