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

#include "sql/zstd_db.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

// ZSTD_c_format / ZSTD_d_format (magicless frames) are in zstd's
// experimental API; the parameters are handled by the shared library too
#define ZSTD_STATIC_LINKING_ONLY
#include <zdict.h>
#include <zstd.h>

#include "common/util/logger.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace zstd_db {

namespace {

const char* const kBlobSuffix = "_blob";
const char* const kCodecTable = "zstd_blob_columns";
const char* const kTablesTable = "zstd_blob_tables";
// sqlite-zstd trains a dictionary of 1% of a column's bytes and none (the
// column stays as is) when that is under min_dict_size_bytes_for_training;
// TopMSV's converter sets 256, zstd's smallest dictionary
const int64_t kMinDictBytes = 256;

[[noreturn]] void fail(const std::string& msg) {
  LOG_ERROR(msg);
  std::exit(EXIT_FAILURE);
}

// "name" and `name`: SQLite's and sqlite-zstd's quoting of identifiers
std::string q(const std::string& name) {
  std::string s = "\"";
  for (char c : name) {
    s += c;
    if (c == '"') s += '"';
  }
  return s + "\"";
}

std::string bq(const std::string& name) {
  std::string s = "`";
  for (char c : name) {
    s += c;
    if (c == '`') s += '`';
  }
  return s + "`";
}

std::string sqlString(const std::string& s) {
  std::string out = "'";
  for (char c : s) {
    out += c;
    if (c == '\'') out += '\'';
  }
  return out + "'";
}

std::string jsonString(const std::string& s) {
  std::string out = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') out += '\\';
    out += c;
  }
  return out + "\"";
}

void removeFile(const std::string& name) {
  std::error_code ec;
  std::filesystem::remove(name, ec);
}

void removeDbFiles(const std::string& name) {
  for (const char* suffix : {"", "-journal", "-wal", "-shm"}) {
    removeFile(name + suffix);
  }
}

int64_t fileSize(const std::string& name) {
  std::error_code ec;
  auto size = std::filesystem::file_size(name, ec);
  return ec ? 0 : static_cast<int64_t>(size);
}

// The source opened read only: as immutable (no -wal / -shm files next to
// it), unless it has a write-ahead log with content, which an immutable open
// would not read
std::string sourceUri(const std::string& name) {
  std::string uri = sql_util::immutableUri(name);
  if (fileSize(name + "-wal") > 0) {
    uri = uri.substr(0, uri.size() - std::string("immutable=1").size()) +
          "mode=ro";
  }
  return uri;
}

std::string columnText(sqlite3_stmt* stmt, int col) {
  const unsigned char* t = sqlite3_column_text(stmt, col);
  return t == nullptr ? "" : reinterpret_cast<const char*>(t);
}

sqlite3* openNew(const std::string& name) {
  removeDbFiles(name);
  sqlite3* db = nullptr;
  if (sqlite3_open_v2(name.c_str(), &db,
                      SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE |
                          SQLITE_OPEN_URI,
                      nullptr) != SQLITE_OK) {
    fail("Cannot open the database " + name + ": " + sqlite3_errmsg(db));
  }
  // a new file under a temporary name: no journal is needed; the rows are
  // copied as they are (TopDIA's MS/MS windows have ms1_id -1)
  sql_util::execSql(db, "PRAGMA journal_mode = OFF;");
  sql_util::execSql(db, "PRAGMA synchronous = OFF;");
  sql_util::execSql(db, "PRAGMA foreign_keys = OFF;");
  return db;
}

void attachSource(sqlite3* db, const std::string& src_name) {
  sqlite3_stmt* stmt = sql_util::prepareSql(db, "ATTACH DATABASE ? AS src;");
  std::string uri = sourceUri(src_name);
  sqlite3_bind_text(stmt, 1, uri.c_str(), -1, SQLITE_TRANSIENT);
  if (sqlite3_step(stmt) != SQLITE_DONE) {
    fail("Cannot open the database " + src_name + ": " + sqlite3_errmsg(db));
  }
  sqlite3_finalize(stmt);
}

// VACUUMs the new file, closes it and renames it to dest_name (replacing an
// earlier file, also the source when compressing in place)
void finish(sqlite3* db, const std::string& tmp_name,
            const std::string& dest_name) {
  sql_util::execSql(db, "PRAGMA journal_mode = DELETE;");
  sql_util::execSql(db, "PRAGMA temp_store = FILE;");
  sql_util::execSql(db, "VACUUM;");
  if (sqlite3_close(db) != SQLITE_OK) {
    fail("Cannot close the database " + tmp_name);
  }
  removeDbFiles(dest_name);
  std::error_code ec;
  std::filesystem::rename(tmp_name, dest_name, ec);
  if (ec) {
    fail("Cannot rename " + tmp_name + " to " + dest_name + ": " +
         ec.message());
  }
}

struct TableDef {
  std::string name;
  std::string sql;
};

std::vector<TableDef> sourceTables(sqlite3* db) {
  std::vector<TableDef> tables;
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db,
      "SELECT name, sql FROM src.sqlite_master WHERE type = 'table' AND name "
      "NOT LIKE 'sqlite_%' ORDER BY rowid;");
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    tables.push_back({columnText(stmt, 0), columnText(stmt, 1)});
  }
  sqlite3_finalize(stmt);
  return tables;
}

bool withoutRowid(const std::string& create_sql) {
  std::string upper = create_sql;
  std::transform(upper.begin(), upper.end(), upper.begin(),
                 [](unsigned char c) { return std::toupper(c); });
  size_t pos = upper.rfind("WITHOUT");
  return pos != std::string::npos &&
         upper.find("ROWID", pos) != std::string::npos;
}

bool declaredInteger(const std::string& type) {
  std::string upper = type;
  std::transform(upper.begin(), upper.end(), upper.begin(),
                 [](unsigned char c) { return std::toupper(c); });
  return upper.find("INT") != std::string::npos;
}

// ------------------------------------------------------------------ codecs

enum class Codec { kF4, kF8, kI1, kI2, kI4, kSeq, kSeqg };

const char* codecName(Codec c) {
  switch (c) {
    case Codec::kF4: return "f4";
    case Codec::kF8: return "f8";
    case Codec::kI1: return "i1";
    case Codec::kI2: return "i2";
    case Codec::kI4: return "i4";
    case Codec::kSeq: return "seq";
    case Codec::kSeqg: return "seqg";
  }
  return "";
}

bool parseCodec(const std::string& s, Codec* c) {
  static const std::vector<std::pair<std::string, Codec>> names = {
      {"f4", Codec::kF4}, {"f8", Codec::kF8},   {"i1", Codec::kI1},
      {"i2", Codec::kI2}, {"i4", Codec::kI4},   {"seq", Codec::kSeq},
      {"seqg", Codec::kSeqg}};
  for (const auto& p : names) {
    if (p.first == s) {
      *c = p.second;
      return true;
    }
  }
  return false;
}

size_t codecWidth(Codec c) {
  switch (c) {
    case Codec::kF4: return 4;
    case Codec::kF8: return 8;
    case Codec::kI1: return 1;
    case Codec::kI2: return 2;
    case Codec::kI4: return 4;
    default: return 0;
  }
}

bool isIntCodec(Codec c) {
  return c != Codec::kF4 && c != Codec::kF8;
}

bool isStored(Codec c) { return c != Codec::kSeq && c != Codec::kSeqg; }

// A column's values over a table, to choose its codec (TopMSV's
// src/server/zstd/codec.ts CodecChooser, the same rules)
class Chooser {
 public:
  void nextSpectrum() {
    pos_ = 0;
    first_ = true;
  }

  void add(double v, bool integral, double group) {
    if (all_int_ && !integral) all_int_ = false;
    if (all_f4_ &&
        static_cast<double>(static_cast<float>(v)) != v) {
      all_f4_ = false;
    }
    if (v < min_) min_ = v;
    if (v > max_) max_ = v;
    if (seq_ && v != static_cast<double>(pos_)) seq_ = false;
    pos_++;
    if (first_ || group != group_) {
      group_ = group;
      gpos_ = 0;
      first_ = false;
    }
    if (seqg_ && v != static_cast<double>(gpos_)) seqg_ = false;
    gpos_++;
  }

  Codec codec(bool declared_integer, bool id_column, bool grouped) const {
    if (declared_integer && all_int_) {
      if (id_column && seq_) return Codec::kSeq;
      if (id_column && grouped && seqg_) return Codec::kSeqg;
      double lo = min_ == kInf ? 0 : min_;
      double hi = max_ == -kInf ? 0 : max_;
      if (lo >= -128 && hi <= 127) return Codec::kI1;
      if (lo >= -32768 && hi <= 32767) return Codec::kI2;
      if (lo >= -2147483648.0 && hi <= 2147483647.0) return Codec::kI4;
      std::ostringstream msg;
      msg << "integer values out of the 32-bit range (" << lo << " - " << hi
          << ")";
      fail(msg.str());
    }
    return all_f4_ ? Codec::kF4 : Codec::kF8;
  }

  double min() const { return min_ == kInf ? 0 : min_; }
  double max() const { return max_ == -kInf ? 0 : max_; }

 private:
  static constexpr double kInf = std::numeric_limits<double>::infinity();
  bool all_int_ = true;
  bool all_f4_ = true;
  double min_ = kInf;
  double max_ = -kInf;
  bool seq_ = true;
  int64_t pos_ = 0;
  bool seqg_ = true;
  int64_t gpos_ = 0;
  bool first_ = true;
  double group_ = 0;
};

// The value of a source cell as a number; exits on a NULL, a text or a blob,
// which the arrays cannot store
double cellValue(sqlite3_stmt* stmt, int col, const std::string& table,
                 const std::string& column, bool* integral) {
  int type = sqlite3_column_type(stmt, col);
  if (type == SQLITE_INTEGER) {
    *integral = true;
    return static_cast<double>(sqlite3_column_int64(stmt, col));
  }
  if (type == SQLITE_FLOAT) {
    double v = sqlite3_column_double(stmt, col);
    *integral = std::isfinite(v) && std::floor(v) == v;
    return v;
  }
  fail("unsupported value in " + table + "." + column +
       " (only non-null numbers are stored)");
}

void appendValue(std::string* bytes, Codec c, double v) {
  // the arrays are little-endian, as the platforms the suite builds on
  switch (c) {
    case Codec::kF4: {
      float f = static_cast<float>(v);
      bytes->append(reinterpret_cast<const char*>(&f), 4);
      break;
    }
    case Codec::kF8:
      bytes->append(reinterpret_cast<const char*>(&v), 8);
      break;
    case Codec::kI1: {
      int8_t i = static_cast<int8_t>(v);
      bytes->append(reinterpret_cast<const char*>(&i), 1);
      break;
    }
    case Codec::kI2: {
      int16_t i = static_cast<int16_t>(v);
      bytes->append(reinterpret_cast<const char*>(&i), 2);
      break;
    }
    case Codec::kI4: {
      int32_t i = static_cast<int32_t>(v);
      bytes->append(reinterpret_cast<const char*>(&i), 4);
      break;
    }
    default:
      break;
  }
}

// ---------------------------------------------------------------- zstd

// one array as sqlite-zstd's zstd_compress_col writes it (src/basic.rs of
// v0.3.5: zstd::bulk::Compressor with the column's prepared dictionary, the
// source size pledged, and with "compact": no checksum, content size,
// dictionary id or magic bytes)
std::string compressValue(ZSTD_CCtx* cctx, const ZSTD_CDict* cdict,
                          const std::string& src) {
  ZSTD_CCtx_reset(cctx, ZSTD_reset_session_and_parameters);
  size_t rc = ZSTD_CCtx_refCDict(cctx, cdict);
  if (!ZSTD_isError(rc)) rc = ZSTD_CCtx_setPledgedSrcSize(cctx, src.size());
  if (!ZSTD_isError(rc)) rc = ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 0);
  if (!ZSTD_isError(rc)) rc = ZSTD_CCtx_setParameter(cctx, ZSTD_c_contentSizeFlag, 0);
  if (!ZSTD_isError(rc)) rc = ZSTD_CCtx_setParameter(cctx, ZSTD_c_dictIDFlag, 0);
  if (!ZSTD_isError(rc)) {
    rc = ZSTD_CCtx_setParameter(cctx, ZSTD_c_format, ZSTD_f_zstd1_magicless);
  }
  if (ZSTD_isError(rc)) {
    fail(std::string("zstd parameters: ") + ZSTD_getErrorName(rc));
  }
  std::string out(ZSTD_compressBound(src.size()), '\0');
  size_t size = ZSTD_compress2(cctx, out.data(), out.size(), src.data(),
                               src.size());
  if (ZSTD_isError(size)) {
    fail(std::string("zstd compression: ") + ZSTD_getErrorName(size));
  }
  out.resize(size);
  return out;
}

// sqlite-zstd's objects for the BLOB table blob_name (the statements the
// extension's zstd_enable_transparent runs, src/transparent.rs of v0.3.5)
void writeTransparentObjects(sqlite3* db, const std::string& blob_name,
                             const std::vector<std::string>& stored) {
  std::string zstd_name = "_" + blob_name + "_zstd";
  for (const std::string& c : stored) {
    sql_util::execSql(db, "CREATE INDEX " + q("_" + blob_name + "_" + c +
                              "_dict_idx") +
                              " ON " + q(zstd_name) + "(" +
                              q("_" + c + "_dict") + ");");
  }
  std::string select = "select `spec_id`, `n`";
  for (const std::string& c : stored) {
    select += ", zstd_decompress_col(" + bq(c) + ", 0, " +
              bq("_" + c + "_dict") + ", true) as " + bq(c);
  }
  sql_util::execSql(db, "\n        create view " + bq(blob_name) +
                            " as\n            " + select +
                            "\n            from " + bq(zstd_name) +
                            "\n        ");
  std::string insert_cols = "spec_id, n";
  std::string insert_vals = "new.`spec_id`,\nnew.`n`";
  for (const std::string& c : stored) {
    insert_cols += ", " + c + ", _" + c + "_dict";
    insert_vals += ",\nnew." + bq(c) + " as " + bq(c) + ", null as " +
                   bq("_" + c + "_dict");
  }
  sql_util::execSql(
      db, "\n            create trigger " + bq(blob_name + "_insert_trigger") +
              "\n            instead of insert on " + bq(blob_name) +
              "\n            for each row\n            begin\n"
              "                insert into " +
              bq(zstd_name) + "(" + insert_cols + ") select " + insert_vals +
              ";\n            end;\n            ");
  sql_util::execSql(
      db, "\n            create trigger " + bq(blob_name + "_delete_trigger") +
              "\n            instead of delete on " + bq(blob_name) +
              "\n            for each row\n            begin\n"
              "                delete from " +
              bq(zstd_name) +
              " where old.`spec_id` = `spec_id`;\n            end;\n"
              "            ");
  std::vector<std::string> all = {"spec_id", "n"};
  all.insert(all.end(), stored.begin(), stored.end());
  for (const std::string& c : all) {
    bool data = c != "spec_id" && c != "n";
    std::string set = bq(c) + " = new." + bq(c);
    if (data) set += ", " + bq("_" + c + "_dict") + " = null";
    sql_util::execSql(
        db, "\n                create trigger " +
                bq(blob_name + "_update_" + c + "_trigger") +
                "\n                instead of update of " + bq(c) + " on " +
                bq(blob_name) +
                "\n                for each row\n                begin\n"
                "                    update " +
                bq(zstd_name) + " set " + set +
                " where old.`spec_id` = `spec_id`;\n                end;\n"
                "                ");
  }
}

struct SpectrumColumn {
  std::string name;
  std::string type;
  int pk;
};

}  // namespace

const std::vector<std::string>& spectrumTables() {
  static const std::vector<std::string> tables = {
      "ms1_peak", "ms2_peak", "ms1_env", "ms1_env_peak", "ms2_env",
      "ms2_env_peak"};
  return tables;
}

bool isSqliteFile(const std::string& name) {
  std::ifstream in(name, std::ios::binary);
  char header[16] = {0};
  if (!in.read(header, 16)) {
    return false;
  }
  return std::memcmp(header, "SQLite format 3\0", 16) == 0;
}

bool isCompressed(sqlite3* db, const std::string& schema) {
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db, "SELECT 1 FROM " + schema +
              ".sqlite_master WHERE type = 'table' AND name = "
              "'_zstd_configs';");
  bool found = sqlite3_step(stmt) == SQLITE_ROW;
  sqlite3_finalize(stmt);
  return found;
}

bool isCompressed(const std::string& name) {
  if (!isSqliteFile(name)) {
    return false;
  }
  sqlite3* db = nullptr;
  if (sqlite3_open_v2(sourceUri(name).c_str(), &db,
                      SQLITE_OPEN_READONLY | SQLITE_OPEN_URI,
                      nullptr) != SQLITE_OK) {
    sqlite3_close(db);
    return false;
  }
  bool found = isCompressed(db, "main");
  sqlite3_close(db);
  return found;
}

Stat compressDb(const std::string& src_name, const std::string& dest_name,
                int level) {
  auto start = std::chrono::steady_clock::now();
  if (level < 1 || level > 22) {
    fail("The zstd level must be 1 - 22.");
  }
  if (!isSqliteFile(src_name)) {
    fail(src_name + " is not an SQLite database.");
  }
  if (isCompressed(src_name)) {
    fail(src_name + " is compressed already.");
  }
  Stat stat;
  stat.bytes_before = fileSize(src_name);
  std::string tmp_name = dest_name + ".tmp";
  sqlite3* db = openNew(tmp_name);
  attachSource(db, src_name);
  std::vector<TableDef> tables = sourceTables(db);
  const std::vector<std::string>& spec_tables = spectrumTables();
  auto isSpectrumTable = [&](const std::string& name) {
    return std::find(spec_tables.begin(), spec_tables.end(), name) !=
           spec_tables.end();
  };

  // 1. the other tables as they are, with their explicit indexes
  sql_util::execSql(db, "BEGIN;");
  std::set<std::string> copied;
  for (const TableDef& t : tables) {
    if (isSpectrumTable(t.name)) continue;
    sql_util::execSql(db, t.sql + ";");
    sql_util::execSql(db, "INSERT INTO main." + q(t.name) +
                              " SELECT * FROM src." + q(t.name) +
                              (withoutRowid(t.sql) ? "" : " ORDER BY rowid") +
                              ";");
    copied.insert(t.name);
  }
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db,
      "SELECT tbl_name, sql FROM src.sqlite_master WHERE type = 'index' AND "
      "sql IS NOT NULL ORDER BY rowid;");
  std::vector<std::string> index_sqls;
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    if (copied.count(columnText(stmt, 0)) > 0) {
      index_sqls.push_back(columnText(stmt, 1));
    }
  }
  sqlite3_finalize(stmt);
  for (const std::string& sql : index_sqls) {
    sql_util::execSql(db, sql + ";");
  }
  sql_util::execSql(
      db, std::string("CREATE TABLE ") + kCodecTable +
              "(tbl TEXT NOT NULL, col TEXT NOT NULL, pos INTEGER NOT NULL, "
              "decl_type TEXT NOT NULL, codec TEXT NOT NULL, group_col TEXT, "
              "min_value REAL, max_value REAL, PRIMARY KEY(tbl, col));");
  sql_util::execSql(db, std::string("CREATE TABLE ") + kTablesTable +
                            "(tbl TEXT PRIMARY KEY, create_sql TEXT NOT NULL);");
  sql_util::execSql(db, "COMMIT;");

  // 2. the per-spectrum tables as BLOB rows (the arrays as they are)
  struct BlobTable {
    std::string name;
    std::vector<std::string> stored;
  };
  std::vector<BlobTable> blob_tables;
  for (const TableDef& t : tables) {
    if (!isSpectrumTable(t.name)) continue;
    std::vector<SpectrumColumn> cols;
    stmt = sql_util::prepareSql(db, "PRAGMA src.table_info(" + q(t.name) + ");");
    while (sqlite3_step(stmt) == SQLITE_ROW) {
      cols.push_back({columnText(stmt, 1), columnText(stmt, 2),
                      sqlite3_column_int(stmt, 5)});
    }
    sqlite3_finalize(stmt);
    int spec_col = -1;
    for (size_t i = 0; i < cols.size(); i++) {
      if (cols[i].name == "spec_id") spec_col = static_cast<int>(i);
    }
    if (spec_col < 0) {
      fail("The table " + t.name + " of " + src_name + " has no spec_id.");
    }
    // the rows in the source's order: rowid order, or for a WITHOUT ROWID
    // table (TopFD's merged ms1_peak, key LEVEL, spec_id, peak_id) the key
    // without LEVEL
    std::string order = "spec_id, rowid";
    if (withoutRowid(t.sql)) {
      std::vector<std::pair<int, std::string>> key;
      for (const SpectrumColumn& c : cols) {
        std::string upper = c.name;
        std::transform(upper.begin(), upper.end(), upper.begin(),
                       [](unsigned char ch) { return std::toupper(ch); });
        if (c.pk > 0 && c.name != "spec_id" && upper != "LEVEL") {
          key.push_back({c.pk, q(c.name)});
        }
      }
      std::sort(key.begin(), key.end());
      order = "spec_id";
      for (const auto& k : key) order += ", " + k.second;
    }
    std::string select = "SELECT ";
    for (size_t i = 0; i < cols.size(); i++) {
      select += (i > 0 ? ", " : "") + q(cols[i].name);
    }
    select += " FROM src." + q(t.name) + " ORDER BY " + order + ";";
    // the data columns, and the group column of an envelope peak's peak_id
    std::vector<size_t> data;
    std::vector<int> group_of;
    int env_col = -1;
    for (size_t i = 0; i < cols.size(); i++) {
      if (cols[i].name == "env_id") env_col = static_cast<int>(i);
    }
    for (size_t i = 0; i < cols.size(); i++) {
      if (static_cast<int>(i) == spec_col) continue;
      data.push_back(i);
      group_of.push_back(cols[i].name == "peak_id" ? env_col : -1);
    }
    // pass 1: the codecs
    std::vector<Chooser> choosers(data.size());
    stmt = sql_util::prepareSql(db, select);
    bool have_spec = false;
    int64_t spec = 0;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
      int64_t s = sqlite3_column_int64(stmt, spec_col);
      if (!have_spec || s != spec) {
        have_spec = true;
        spec = s;
        for (Chooser& c : choosers) c.nextSpectrum();
      }
      for (size_t k = 0; k < data.size(); k++) {
        bool integral = false;
        double v = cellValue(stmt, static_cast<int>(data[k]), t.name,
                             cols[data[k]].name, &integral);
        double group = 0;
        if (group_of[k] >= 0) {
          bool g_int = false;
          group = cellValue(stmt, group_of[k], t.name, "env_id", &g_int);
        }
        choosers[k].add(v, integral, group);
      }
    }
    sqlite3_finalize(stmt);
    std::vector<Codec> codecs(data.size());
    BlobTable bt;
    bt.name = t.name + kBlobSuffix;
    for (size_t k = 0; k < data.size(); k++) {
      const SpectrumColumn& c = cols[data[k]];
      codecs[k] = choosers[k].codec(declaredInteger(c.type),
                                    c.name == "peak_id" || c.name == "env_id",
                                    group_of[k] >= 0);
      if (isStored(codecs[k])) bt.stored.push_back(c.name);
    }
    // the table as sqlite-zstd leaves it: renamed _<blob>_zstd, a dictionary
    // id column per array column
    std::string create = "CREATE TABLE " + q("_" + bt.name + "_zstd") +
                         "(spec_id INTEGER PRIMARY KEY, n INTEGER NOT NULL";
    for (const std::string& c : bt.stored) create += ", " + q(c) + " BLOB NOT NULL";
    for (const std::string& c : bt.stored) {
      create += ", " + bq("_" + c + "_dict") +
                " integer default null references _zstd_dicts(id)";
    }
    create += ")";
    sql_util::execSql(db, "BEGIN;");
    sql_util::execSql(db, create + ";");
    sqlite3_stmt* codec_stmt = sql_util::prepareSql(
        db, std::string("INSERT INTO ") + kCodecTable +
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?);");
    for (size_t i = 0; i < cols.size(); i++) {
      sqlite3_bind_text(codec_stmt, 1, t.name.c_str(), -1, SQLITE_TRANSIENT);
      sqlite3_bind_text(codec_stmt, 2, cols[i].name.c_str(), -1,
                        SQLITE_TRANSIENT);
      sqlite3_bind_int(codec_stmt, 3, static_cast<int>(i));
      sqlite3_bind_text(codec_stmt, 4, cols[i].type.c_str(), -1,
                        SQLITE_TRANSIENT);
      auto it = std::find(data.begin(), data.end(), i);
      if (it == data.end()) {
        sqlite3_bind_text(codec_stmt, 5, "key", -1, SQLITE_STATIC);
        sqlite3_bind_null(codec_stmt, 6);
        sqlite3_bind_null(codec_stmt, 7);
        sqlite3_bind_null(codec_stmt, 8);
      } else {
        size_t k = static_cast<size_t>(it - data.begin());
        sqlite3_bind_text(codec_stmt, 5, codecName(codecs[k]), -1,
                          SQLITE_STATIC);
        if (codecs[k] == Codec::kSeqg) {
          sqlite3_bind_text(codec_stmt, 6, "env_id", -1, SQLITE_STATIC);
        } else {
          sqlite3_bind_null(codec_stmt, 6);
        }
        sqlite3_bind_double(codec_stmt, 7, choosers[k].min());
        sqlite3_bind_double(codec_stmt, 8, choosers[k].max());
      }
      sql_util::stepAndReset(db, codec_stmt);
    }
    sqlite3_finalize(codec_stmt);
    sqlite3_stmt* tbl_stmt = sql_util::prepareSql(
        db, std::string("INSERT INTO ") + kTablesTable + " VALUES (?, ?);");
    sqlite3_bind_text(tbl_stmt, 1, t.name.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(tbl_stmt, 2, t.sql.c_str(), -1, SQLITE_TRANSIENT);
    sql_util::stepAndReset(db, tbl_stmt);
    sqlite3_finalize(tbl_stmt);
    // pass 2: one row per spectrum
    std::string ins_sql = "INSERT INTO " + q("_" + bt.name + "_zstd") +
                          "(spec_id, n";
    for (const std::string& c : bt.stored) ins_sql += ", " + q(c);
    ins_sql += ") VALUES (?, ?";
    for (size_t k = 0; k < bt.stored.size(); k++) ins_sql += ", ?";
    ins_sql += ");";
    sqlite3_stmt* ins = sql_util::prepareSql(db, ins_sql);
    std::vector<std::string> arrays(data.size());
    int64_t n = 0;
    auto flush = [&]() {
      if (n == 0) return;
      sqlite3_bind_int64(ins, 1, spec);
      sqlite3_bind_int64(ins, 2, n);
      int p = 3;
      for (size_t k = 0; k < data.size(); k++) {
        if (!isStored(codecs[k])) continue;
        sqlite3_bind_blob(ins, p++, arrays[k].data(),
                          static_cast<int>(arrays[k].size()), SQLITE_STATIC);
      }
      sql_util::stepAndReset(db, ins);
    };
    stmt = sql_util::prepareSql(db, select);
    have_spec = false;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
      int64_t s = sqlite3_column_int64(stmt, spec_col);
      if (!have_spec || s != spec) {
        flush();
        have_spec = true;
        spec = s;
        for (std::string& a : arrays) a.clear();
        n = 0;
      }
      for (size_t k = 0; k < data.size(); k++) {
        if (!isStored(codecs[k])) continue;
        bool integral = false;
        appendValue(&arrays[k], codecs[k],
                    cellValue(stmt, static_cast<int>(data[k]), t.name,
                              cols[data[k]].name, &integral));
      }
      n++;
    }
    flush();
    sqlite3_finalize(stmt);
    sqlite3_finalize(ins);
    sql_util::execSql(db, "COMMIT;");
    blob_tables.push_back(bt);
  }
  sql_util::execSql(db, "DETACH DATABASE src;");

  // 3. sqlite-zstd: a dictionary per BLOB column (trained on all its
  // arrays) and every array compressed with it
  sql_util::execSql(db, "BEGIN;");
  sql_util::execSql(db,
                    "\n        create table if not exists _zstd_dicts (\n"
                    "            id integer primary key autoincrement,\n"
                    "            chooser_key text unique,\n"
                    "            dict blob not null\n        );");
  sql_util::execSql(db,
                    "insert or ignore into _zstd_dicts values (-1, "
                    "'[nodict]', '');");
  sql_util::execSql(db,
                    "\n            create table if not exists _zstd_configs (\n"
                    "                id integer primary key autoincrement,\n"
                    "                config json not null\n            );");
  sqlite3_stmt* cfg_stmt = sql_util::prepareSql(
      db, "INSERT INTO _zstd_configs (config) VALUES (?);");
  for (const BlobTable& bt : blob_tables) {
    for (const std::string& c : bt.stored) {
      std::ostringstream cfg;
      cfg << "{\"table\":" << jsonString(bt.name)
          << ",\"column\":" << jsonString(c)
          << ",\"compression_level\":" << level
          << ",\"dict_chooser\":" << jsonString(sqlString(bt.name + "." + c))
          << ",\"min_dict_size_bytes_for_training\":" << kMinDictBytes << "}";
      std::string text = cfg.str();
      sqlite3_bind_text(cfg_stmt, 1, text.c_str(), -1, SQLITE_TRANSIENT);
      sql_util::stepAndReset(db, cfg_stmt);
    }
  }
  sqlite3_finalize(cfg_stmt);
  sql_util::execSql(db, "COMMIT;");
  ZSTD_CCtx* cctx = ZSTD_createCCtx();
  for (const BlobTable& bt : blob_tables) {
    std::string zstd_name = "_" + bt.name + "_zstd";
    for (const std::string& c : bt.stored) {
      // the column's arrays in rowid (spec_id) order
      std::vector<int64_t> ids;
      std::string samples;
      std::vector<size_t> sizes;
      stmt = sql_util::prepareSql(db, "SELECT spec_id, " + q(c) + " FROM " +
                                          q(zstd_name) + " ORDER BY spec_id;");
      while (sqlite3_step(stmt) == SQLITE_ROW) {
        ids.push_back(sqlite3_column_int64(stmt, 0));
        const char* b = static_cast<const char*>(sqlite3_column_blob(stmt, 1));
        int len = sqlite3_column_bytes(stmt, 1);
        samples.append(b == nullptr ? "" : std::string(b, len));
        sizes.push_back(static_cast<size_t>(len));
      }
      sqlite3_finalize(stmt);
      if (ids.empty()) continue;
      // sqlite-zstd: (total bytes as f32 * dict_size_ratio 0.01) as i64
      int64_t dict_size = static_cast<int64_t>(
          static_cast<float>(samples.size()) * 0.01f);
      if (dict_size < kMinDictBytes) {
        continue;  // stored as is, dictionary id NULL
      }
      std::string dict(static_cast<size_t>(dict_size), '\0');
      size_t got = ZDICT_trainFromBuffer(dict.data(), dict.size(),
                                         samples.data(), sizes.data(),
                                         static_cast<unsigned>(sizes.size()));
      if (ZDICT_isError(got)) {
        LOG_WARN("No zstd dictionary for " << bt.name << "." << c << ": "
                                           << ZDICT_getErrorName(got)
                                           << "; its arrays are stored as is.");
        continue;
      }
      dict.resize(got);
      sql_util::execSql(db, "BEGIN;");
      sqlite3_stmt* dict_stmt = sql_util::prepareSql(
          db, "insert into _zstd_dicts (chooser_key,dict) values (?, ?);");
      std::string key = bt.name + "." + c;
      sqlite3_bind_text(dict_stmt, 1, key.c_str(), -1, SQLITE_TRANSIENT);
      sqlite3_bind_blob(dict_stmt, 2, dict.data(), static_cast<int>(dict.size()),
                        SQLITE_TRANSIENT);
      sql_util::stepAndReset(db, dict_stmt);
      sqlite3_finalize(dict_stmt);
      int64_t dict_id = sqlite3_last_insert_rowid(db);
      ZSTD_CDict* cdict = ZSTD_createCDict(dict.data(), dict.size(), level);
      sqlite3_stmt* upd = sql_util::prepareSql(
          db, "UPDATE " + q(zstd_name) + " SET " + q(c) + " = ?, " +
                  q("_" + c + "_dict") + " = ? WHERE spec_id = ?;");
      size_t offset = 0;
      for (size_t i = 0; i < ids.size(); i++) {
        std::string value = samples.substr(offset, sizes[i]);
        offset += sizes[i];
        std::string packed = compressValue(cctx, cdict, value);
        sqlite3_bind_blob(upd, 1, packed.data(), static_cast<int>(packed.size()),
                          SQLITE_TRANSIENT);
        sqlite3_bind_int64(upd, 2, dict_id);
        sqlite3_bind_int64(upd, 3, ids[i]);
        sql_util::stepAndReset(db, upd);
      }
      sqlite3_finalize(upd);
      ZSTD_freeCDict(cdict);
      sql_util::execSql(db, "COMMIT;");
    }
    sql_util::execSql(db, "BEGIN;");
    writeTransparentObjects(db, bt.name, bt.stored);
    sql_util::execSql(db, "COMMIT;");
  }
  ZSTD_freeCCtx(cctx);
  finish(db, tmp_name, dest_name);
  stat.bytes_after = fileSize(dest_name);
  stat.seconds = std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - start)
                     .count();
  return stat;
}

Stat decompressDb(const std::string& src_name, const std::string& dest_name) {
  auto start = std::chrono::steady_clock::now();
  if (!isSqliteFile(src_name)) {
    fail(src_name + " is not an SQLite database.");
  }
  if (!isCompressed(src_name)) {
    fail(src_name + " is not compressed.");
  }
  Stat stat;
  stat.bytes_before = fileSize(src_name);
  std::string tmp_name = dest_name + ".tmp";
  sqlite3* db = openNew(tmp_name);
  attachSource(db, src_name);
  std::vector<TableDef> tables = sourceTables(db);
  // the per-spectrum tables, in the order they were compressed
  std::vector<std::string> blob_tables;
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db, std::string("SELECT tbl FROM src.") + kCodecTable +
              " GROUP BY tbl ORDER BY MIN(rowid);");
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    blob_tables.push_back(columnText(stmt, 0));
  }
  sqlite3_finalize(stmt);
  std::set<std::string> skip = {kCodecTable, kTablesTable, "_zstd_configs",
                                "_zstd_dicts"};
  for (const std::string& t : blob_tables) {
    skip.insert("_" + t + kBlobSuffix + "_zstd");
  }
  sql_util::execSql(db, "BEGIN;");
  std::set<std::string> copied;
  for (const TableDef& t : tables) {
    if (skip.count(t.name) > 0) continue;
    sql_util::execSql(db, t.sql + ";");
    sql_util::execSql(db, "INSERT INTO main." + q(t.name) +
                              " SELECT * FROM src." + q(t.name) +
                              (withoutRowid(t.sql) ? "" : " ORDER BY rowid") +
                              ";");
    copied.insert(t.name);
  }
  stmt = sql_util::prepareSql(
      db,
      "SELECT tbl_name, sql FROM src.sqlite_master WHERE type = 'index' AND "
      "sql IS NOT NULL ORDER BY rowid;");
  std::vector<std::string> index_sqls;
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    if (copied.count(columnText(stmt, 0)) > 0) {
      index_sqls.push_back(columnText(stmt, 1));
    }
  }
  sqlite3_finalize(stmt);
  for (const std::string& sql : index_sqls) {
    sql_util::execSql(db, sql + ";");
  }
  for (const std::string& t : blob_tables) {
    BlobTableReader reader(db, "src", t);
    sql_util::execSql(db, reader.createSql() + ";");
    const std::vector<std::string>& names = reader.columnNames();
    std::string sql = "INSERT INTO main." + q(t) + " VALUES (";
    for (size_t i = 0; i < names.size(); i++) sql += (i > 0 ? ", ?" : "?");
    sql += ");";
    sqlite3_stmt* ins = sql_util::prepareSql(db, sql);
    while (reader.next()) {
      for (int r = 0; r < reader.rowNum(); r++) {
        for (size_t c = 0; c < names.size(); c++) {
          reader.bind(ins, static_cast<int>(c) + 1, c, r);
        }
        sql_util::stepAndReset(db, ins);
      }
    }
    sqlite3_finalize(ins);
  }
  sql_util::execSql(db, "COMMIT;");
  sql_util::execSql(db, "DETACH DATABASE src;");
  finish(db, tmp_name, dest_name);
  stat.bytes_after = fileSize(dest_name);
  stat.seconds = std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - start)
                     .count();
  return stat;
}

void compressInPlace(const std::string& name, int level) {
  if (!std::filesystem::exists(name)) {
    return;
  }
  std::string file = std::filesystem::path(name).filename().string();
  if (isCompressed(name)) {
    std::cout << "The SQLite database " << file << " is compressed already."
              << std::endl;
    return;
  }
  std::cout << "Compressing the SQLite database " << file << " started."
            << std::endl;
  Stat stat = compressDb(name, name, level);
  std::cout << "Compressing the SQLite database " << file << " finished ("
            << std::fixed << std::setprecision(1) << stat.bytes_before / 1e6
            << " MB -> " << stat.bytes_after / 1e6 << " MB, " << stat.seconds
            << " s)." << std::defaultfloat << std::endl;
}

// ------------------------------------------------------- BlobTableReader

struct BlobTableReader::Column {
  std::string name;
  Codec codec = Codec::kF8;
  int group = -1;            // seqg: the index of the group column
  int data_col = -1;         // the array's column in the row
  int dict_col = -1;         // its dictionary id's column
  std::vector<int64_t> ints;
  std::vector<double> reals;
};

BlobTableReader::BlobTableReader(sqlite3* db, const std::string& schema,
                                 const std::string& table)
    : db_(db), schema_(schema) {
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db, "SELECT name FROM " + schema +
              ".sqlite_master WHERE type = 'table' AND name IN ('" +
              kCodecTable + "', '" + kTablesTable + "');");
  bool has_codecs = false;
  bool has_tables = false;
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    std::string name = columnText(stmt, 0);
    if (name == kCodecTable) has_codecs = true;
    if (name == kTablesTable) has_tables = true;
  }
  sqlite3_finalize(stmt);
  if (!has_codecs) {
    return;
  }
  std::vector<std::string> types;
  std::vector<std::string> group_names;
  stmt = sql_util::prepareSql(
      db, "SELECT col, decl_type, codec, group_col FROM " + schema + "." +
              kCodecTable + " WHERE tbl = ? ORDER BY pos;");
  sqlite3_bind_text(stmt, 1, table.c_str(), -1, SQLITE_TRANSIENT);
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    auto col = std::make_unique<Column>();
    col->name = columnText(stmt, 0);
    types.push_back(columnText(stmt, 1));
    std::string codec = columnText(stmt, 2);
    group_names.push_back(columnText(stmt, 3));
    bool key = codec == "key";
    if (!key && !parseCodec(codec, &col->codec)) {
      fail("unknown codec " + codec + " of " + table + "." + col->name);
    }
    names_.push_back(col->name);
    is_int_.push_back(key || isIntCodec(col->codec));
    if (key) col->codec = Codec::kI4;  // spec_id: the row's key
    col->group = key ? -2 : -1;        // -2 marks spec_id
    cols_.push_back(std::move(col));
  }
  sqlite3_finalize(stmt);
  if (cols_.empty()) {
    return;
  }
  exists_ = true;
  for (size_t i = 0; i < cols_.size(); i++) {
    if (cols_[i]->codec == Codec::kSeqg && cols_[i]->group != -2) {
      for (size_t j = 0; j < names_.size(); j++) {
        if (names_[j] == group_names[i]) cols_[i]->group = static_cast<int>(j);
      }
      if (cols_[i]->group < 0) {
        fail("no group column " + group_names[i] + " for " + table + "." +
             cols_[i]->name);
      }
    }
  }
  if (has_tables) {
    stmt = sql_util::prepareSql(db, "SELECT create_sql FROM " + schema + "." +
                                        kTablesTable + " WHERE tbl = ?;");
    sqlite3_bind_text(stmt, 1, table.c_str(), -1, SQLITE_TRANSIENT);
    if (sqlite3_step(stmt) == SQLITE_ROW) create_sql_ = columnText(stmt, 0);
    sqlite3_finalize(stmt);
  }
  if (create_sql_.empty()) {
    // a file of TopMSV's first converter (toppic_test experiment 171)
    create_sql_ = "CREATE TABLE " + table + "(";
    for (size_t i = 0; i < names_.size(); i++) {
      create_sql_ += (i > 0 ? "," : "") + names_[i] + " " + types[i];
    }
    create_sql_ += ")";
  }
  std::string from = schema + "." + q("_" + table + kBlobSuffix + "_zstd");
  all_stmt_ = sql_util::prepareSql(db, "SELECT * FROM " + from +
                                           " ORDER BY spec_id;");
  one_stmt_ = sql_util::prepareSql(db, "SELECT * FROM " + from +
                                           " WHERE spec_id = ?;");
  // the columns of the stored rows
  int count = sqlite3_column_count(all_stmt_);
  for (int i = 0; i < count; i++) {
    std::string name = sqlite3_column_name(all_stmt_, i);
    for (auto& col : cols_) {
      if (col->name == name) col->data_col = i;
      if ("_" + col->name + "_dict" == name) col->dict_col = i;
    }
  }
  ZSTD_DCtx* dctx = ZSTD_createDCtx();
  ZSTD_DCtx_setParameter(dctx, ZSTD_d_format, ZSTD_f_zstd1_magicless);
  dctx_ = dctx;
}

BlobTableReader::~BlobTableReader() {
  if (all_stmt_ != nullptr) sqlite3_finalize(all_stmt_);
  if (one_stmt_ != nullptr) sqlite3_finalize(one_stmt_);
  for (auto& entry : ddicts_) {
    ZSTD_freeDDict(static_cast<ZSTD_DDict*>(entry.second));
  }
  if (dctx_ != nullptr) ZSTD_freeDCtx(static_cast<ZSTD_DCtx*>(dctx_));
}

const std::string& BlobTableReader::dict(int64_t id) {
  auto it = dicts_.find(id);
  if (it != dicts_.end()) return it->second;
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db_, "SELECT dict FROM " + schema_ + "._zstd_dicts WHERE id = ?;");
  sqlite3_bind_int64(stmt, 1, id);
  if (sqlite3_step(stmt) != SQLITE_ROW) {
    fail("no zstd dictionary " + std::to_string(id) + " in the database");
  }
  const char* b = static_cast<const char*>(sqlite3_column_blob(stmt, 0));
  std::string d = b == nullptr ? "" : std::string(b, sqlite3_column_bytes(stmt, 0));
  sqlite3_finalize(stmt);
  ddicts_[id] = ZSTD_createDDict(d.data(), d.size());
  return dicts_[id] = d;
}

bool BlobTableReader::decodeRow(sqlite3_stmt* stmt) {
  spec_id_ = sqlite3_column_int(stmt, 0);
  for (int i = 0; i < sqlite3_column_count(stmt); i++) {
    if (std::string(sqlite3_column_name(stmt, i)) == "n") {
      n_ = sqlite3_column_int(stmt, i);
    } else if (std::string(sqlite3_column_name(stmt, i)) == "spec_id") {
      spec_id_ = sqlite3_column_int(stmt, i);
    }
  }
  ZSTD_DCtx* dctx = static_cast<ZSTD_DCtx*>(dctx_);
  // the stored arrays first, then the seqg columns, which follow them
  for (int pass = 0; pass < 2; pass++) {
    for (auto& col : cols_) {
      if (col->group == -2) continue;
      bool seqg = col->codec == Codec::kSeqg;
      if (seqg != (pass == 1)) continue;
      col->ints.clear();
      col->reals.clear();
      if (col->codec == Codec::kSeq) {
        for (int i = 0; i < n_; i++) col->ints.push_back(i);
        continue;
      }
      if (seqg) {
        const std::vector<int64_t>& g = cols_[col->group]->ints;
        int64_t k = 0;
        for (int i = 0; i < n_; i++) {
          k = (i > 0 && g[i] == g[i - 1]) ? k + 1 : 0;
          col->ints.push_back(k);
        }
        continue;
      }
      const char* src =
          static_cast<const char*>(sqlite3_column_blob(stmt, col->data_col));
      size_t src_len = static_cast<size_t>(sqlite3_column_bytes(stmt, col->data_col));
      size_t width = codecWidth(col->codec);
      std::string raw;
      if (col->dict_col < 0 ||
          sqlite3_column_type(stmt, col->dict_col) == SQLITE_NULL) {
        raw.assign(src == nullptr ? "" : std::string(src, src_len));
      } else {
        int64_t id = sqlite3_column_int64(stmt, col->dict_col);
        const ZSTD_DDict* ddict = nullptr;
        if (id != -1) {
          dict(id);
          ddict = static_cast<const ZSTD_DDict*>(ddicts_[id]);
        }
        ZSTD_DCtx_refDDict(dctx, ddict);
        raw.resize(static_cast<size_t>(n_) * width);
        size_t got = ZSTD_decompressDCtx(dctx, raw.data(), raw.size(), src,
                                         src_len);
        if (ZSTD_isError(got) || got != raw.size()) {
          fail("cannot decompress " + col->name + " of spectrum " +
               std::to_string(spec_id_) +
               (ZSTD_isError(got) ? std::string(": ") + ZSTD_getErrorName(got)
                                  : std::string(": wrong size")));
        }
      }
      if (raw.size() != static_cast<size_t>(n_) * width) {
        fail("wrong array size of " + col->name + " of spectrum " +
             std::to_string(spec_id_));
      }
      for (int i = 0; i < n_; i++) {
        const char* p = raw.data() + static_cast<size_t>(i) * width;
        switch (col->codec) {
          case Codec::kF4: {
            float f;
            std::memcpy(&f, p, 4);
            col->reals.push_back(static_cast<double>(f));
            break;
          }
          case Codec::kF8: {
            double d;
            std::memcpy(&d, p, 8);
            col->reals.push_back(d);
            break;
          }
          case Codec::kI1: {
            int8_t v;
            std::memcpy(&v, p, 1);
            col->ints.push_back(v);
            break;
          }
          case Codec::kI2: {
            int16_t v;
            std::memcpy(&v, p, 2);
            col->ints.push_back(v);
            break;
          }
          case Codec::kI4: {
            int32_t v;
            std::memcpy(&v, p, 4);
            col->ints.push_back(v);
            break;
          }
          default:
            break;
        }
      }
    }
  }
  return true;
}

bool BlobTableReader::next() {
  if (!exists_ || sqlite3_step(all_stmt_) != SQLITE_ROW) {
    return false;
  }
  return decodeRow(all_stmt_);
}

bool BlobTableReader::read(int spec_id) {
  if (!exists_) return false;
  sqlite3_reset(one_stmt_);
  sqlite3_bind_int(one_stmt_, 1, spec_id);
  bool found = sqlite3_step(one_stmt_) == SQLITE_ROW && decodeRow(one_stmt_);
  sqlite3_reset(one_stmt_);
  return found;
}

int64_t BlobTableReader::intValue(size_t col, int row) const {
  const Column& c = *cols_[col];
  if (c.group == -2) return spec_id_;
  if (!c.ints.empty()) return c.ints[row];
  return static_cast<int64_t>(c.reals[row]);
}

double BlobTableReader::realValue(size_t col, int row) const {
  const Column& c = *cols_[col];
  if (c.group == -2) return spec_id_;
  if (!c.reals.empty()) return c.reals[row];
  return static_cast<double>(c.ints[row]);
}

void BlobTableReader::bind(sqlite3_stmt* stmt, int idx, size_t col,
                           int row) const {
  if (is_int_[col]) {
    sqlite3_bind_int64(stmt, idx, intValue(col, row));
  } else {
    sqlite3_bind_double(stmt, idx, realValue(col, row));
  }
}

}  // namespace zstd_db

}  // namespace toppic
