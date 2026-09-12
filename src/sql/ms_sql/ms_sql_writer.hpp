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

#ifndef TOPPIC_SQL_MS_SQL_MS_SQL_WRITER_HPP_
#define TOPPIC_SQL_MS_SQL_MS_SQL_WRITER_HPP_

#include <sqlite3.h>

#include <string>
#include <vector>

namespace toppic {

// One MS1 peak of an LC-MS map.
struct MsSqlPeak {
  double mz = 0.0;
  double inte = 0.0;
  int rt = 0;     // retention time in milliseconds
  int color = 0;  // intensity bucket, 0 .. MsSqlWriter::COLOR_NUM - 1
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

// Writes the MS1 peaks of an LC-MS map into the SQLite database read by the
// 3D visualization.
//
// Layout: table PEAKS0 holds every MS1 peak, and PEAKS1, PEAKS2, ... hold
// progressively down-sampled copies of it. The map is cut into a grid of
// rt_size x mz_size blocks and only the most intense peak of each block is
// kept; both block sizes double for each further layer, until a layer has
// fewer than MIN_LAYER_PEAKS peaks. Row i of the CONFIG table is the range of
// PEAKSi, and each layer is indexed on (RETENTIONTIME, MZ). Every peak carries
// a COLOR bucket assigned from its intensity on a log scale over the whole
// map.
class MsSqlWriter {
 public:
  // Opens the database file, replacing an existing one.
  explicit MsSqlWriter(const std::string& db_file_name);

  ~MsSqlWriter();

  MsSqlWriter(const MsSqlWriter&) = delete;
  MsSqlWriter& operator=(const MsSqlWriter&) = delete;

  // Writes all the tables in one transaction. mz_size is the m/z width of a
  // grid block of the first down-sampled layer (capped at the map's m/z
  // range); the block's retention-time height is the mean MS1 scan interval
  // divided by rt_divider.
  void write(std::vector<MsSqlPeak> peaks, int ms1_scan_num, double mz_size,
             double rt_divider);

  static constexpr int COLOR_NUM = 7;
  static constexpr double DEFAULT_MZ_SIZE = 0.05;
  static constexpr double DEFAULT_RT_DIVIDER = 1.0;

 private:
  void createConfigTable();
  void insertConfig(const MsSqlRange& range);
  void createLayerTable(int layer);
  void insertLayerPeaks(const std::vector<MsSqlPeak>& peaks, int layer);
  void createLayerIndex(int layer);

  sqlite3* db_ = nullptr;

  // A layer with fewer peaks than this is the last one.
  static constexpr int MIN_LAYER_PEAKS = 3000;
  // Growth of the grid blocks from one layer to the next.
  static constexpr double MZ_SCALE = 2.0;
  static constexpr double RT_SCALE = 2.0;
};

}  // namespace toppic

#endif
