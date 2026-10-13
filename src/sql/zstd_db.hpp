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

#ifndef TOPPIC_SQL_ZSTD_DB_HPP_
#define TOPPIC_SQL_ZSTD_DB_HPP_

#include <sqlite3.h>

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace toppic {

namespace zstd_db {

// Compressed SQLite databases (toppic_test experiment 172): the format of
// TopMSV's converter (npm run compress, toppic_test experiment 171), which
// the sqlite-zstd extension (github.com/phiresky/sqlite-zstd, v0.3.5) reads,
// written here with libzstd alone (the suite's SQLite cannot load the
// extension, which needs SQLite >= 3.49.1).
//
// Each per-spectrum table T (ms1_peak, ms2_peak, ms1_env, ms1_env_peak,
// ms2_env, ms2_env_peak) is one row per spectrum, (spec_id, n, an array per
// column) in the source's row order: float32 when every value of the column
// is exactly a float32, else float64; integers in the smallest type; an id
// column left out when it is 0, 1, ... in every spectrum (seq), or in every
// envelope for an envelope peak's peak_id (seqg). zstd_blob_columns has each
// column's codec, position, declared type, group column and range;
// zstd_blob_tables the source's CREATE statement. Every other table is
// copied as it is. The arrays are stored as sqlite-zstd stores them: table
// _T_blob_zstd with a column _<c>_dict per array column (the id of its
// dictionary in _zstd_dicts, NULL for an array stored as is), behind the view
// T_blob that the extension decompresses, with its triggers, _zstd_configs
// and the dictionary indexes. Each column has its own zstd dictionary
// trained on all its arrays (1% of their bytes; none, and the arrays stored
// as is, when that is under 256 bytes); each array is one "compact" zstd
// frame (no magic bytes, dictionary id, content size or checksum) at the
// given level.

// The per-spectrum tables.
const std::vector<std::string>& spectrumTables();

// Whether the database file name is compressed (has sqlite-zstd's table
// _zstd_configs); opened read only as immutable. False for a file that does
// not exist or is not an SQLite database.
bool isCompressed(const std::string& name);

// Whether the database schema ("main", or an attached one) of db is
// compressed.
bool isCompressed(sqlite3* db, const std::string& schema = "main");

// Whether the file is an SQLite database (its header).
bool isSqliteFile(const std::string& name);

struct Stat {
  int64_t bytes_before = 0;
  int64_t bytes_after = 0;
  double seconds = 0;
};

// Writes dest_name, the compressed copy of src_name, at the zstd level (1 -
// 22). dest_name may be src_name (compression in place). The copy is written
// under dest_name + ".tmp" and renamed, so that an interrupted run leaves
// no partial file; src_name is opened read only (as immutable, unless it has
// a write-ahead log with content). Exits with an error message when
// src_name is compressed already or cannot be read, or dest_name cannot be
// written.
Stat compressDb(const std::string& src_name, const std::string& dest_name,
                int level = 3);

// Writes dest_name (may be src_name), the database src_name with its
// per-spectrum tables back as rows: every table as it was before
// compression (CREATE statement, rows in order, values and types), without
// sqlite-zstd's objects. Exits with an error message when src_name is not
// compressed or cannot be read, or dest_name cannot be written.
Stat decompressDb(const std::string& src_name, const std::string& dest_name);

// Compresses the database name in place, with lines "Compressing the SQLite
// database ... started." and "... finished (a MB -> b MB, t s)."; a
// database already compressed is left as it is ("... is compressed
// already."). The default of topfd, topdia, toppic and topmg (their option
// --disable-zstd skips it).
void compressInPlace(const std::string& name, int level = 3);

// One per-spectrum table of a compressed database, decoded spectrum by
// spectrum. Values come as int64 (integer codecs) or double (float codecs),
// the types the source held (zstd_blob_columns' declared types).
class BlobTableReader {
 public:
  // The table of the database schema of db; exists() is false when the
  // schema does not store it compressed.
  BlobTableReader(sqlite3* db, const std::string& schema,
                  const std::string& table);
  ~BlobTableReader();
  BlobTableReader(const BlobTableReader&) = delete;
  BlobTableReader& operator=(const BlobTableReader&) = delete;

  bool exists() const { return exists_; }

  // The source table's columns in order (spec_id included), whether each was
  // stored with an integer codec, and its CREATE statement (rebuilt from the
  // names and declared types for a file without zstd_blob_tables).
  const std::vector<std::string>& columnNames() const { return names_; }
  bool isInteger(size_t col) const { return is_int_[col]; }
  const std::string& createSql() const { return create_sql_; }

  // The next spectrum in spec_id order; false after the last.
  bool next();
  // The spectrum spec_id; false when it has no rows.
  bool read(int spec_id);

  int specId() const { return spec_id_; }
  int rowNum() const { return n_; }
  int64_t intValue(size_t col, int row) const;
  double realValue(size_t col, int row) const;
  // Binds row's value of col to the statement's parameter idx.
  void bind(sqlite3_stmt* stmt, int idx, size_t col, int row) const;

 private:
  struct Column;
  bool decodeRow(sqlite3_stmt* stmt);
  const std::string& dict(int64_t id);

  sqlite3* db_;
  std::string schema_;
  bool exists_ = false;
  std::vector<std::string> names_;
  std::vector<bool> is_int_;
  std::string create_sql_;
  std::vector<std::unique_ptr<Column>> cols_;
  sqlite3_stmt* all_stmt_ = nullptr;
  sqlite3_stmt* one_stmt_ = nullptr;
  std::map<int64_t, std::string> dicts_;
  void* dctx_ = nullptr;
  std::map<int64_t, void*> ddicts_;
  int spec_id_ = -1;
  int n_ = 0;
};

}  // namespace zstd_db

}  // namespace toppic

#endif  // TOPPIC_SQL_ZSTD_DB_HPP_
