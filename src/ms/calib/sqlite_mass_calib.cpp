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

#include "ms/calib/sqlite_mass_calib.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

#include <sqlite3.h>

#include "common/util/logger.hpp"
#include "common/util/time_util.hpp"
#include "common/util/version.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace sqlite_mass_calib {

namespace {

// toppic's own tables, which every run writes again (prsm_sql_output)
const std::vector<std::string> kToppicTables = {
    "prsm",          "proteoform",         "protein",
    "prsm_mass_shift", "prsm_protein_match", "fasta_seq"};

// The user data of the SQL functions below: the model, the spectrum level and
// the retention time (minutes) of every spectrum id of that level.
struct LevelContext {
  const MassCalibModel* model = nullptr;
  int ms_level = 1;
  std::unordered_map<int, double> rt_mins;
};

bool findRt(sqlite3_context* ctx, const LevelContext* c, int spec_id,
            double* rt_min) {
  auto it = c->rt_mins.find(spec_id);
  if (it == c->rt_mins.end()) {
    sqlite3_result_error(ctx, "a spectrum id without a spectrum", -1);
    return false;
  }
  *rt_min = it->second;
  return true;
}

// calib_mz(spec_id, mz)
void calibMz(sqlite3_context* ctx, int, sqlite3_value** argv) {
  const LevelContext* c =
      static_cast<const LevelContext*>(sqlite3_user_data(ctx));
  double rt_min = 0.0;
  if (!findRt(ctx, c, sqlite3_value_int(argv[0]), &rt_min)) {
    return;
  }
  sqlite3_result_double(ctx, c->model->correctMz(
                                 c->ms_level, sqlite3_value_double(argv[1]),
                                 rt_min));
}

// calib_mass(spec_id, mass, charge)
void calibMass(sqlite3_context* ctx, int, sqlite3_value** argv) {
  const LevelContext* c =
      static_cast<const LevelContext*>(sqlite3_user_data(ctx));
  double rt_min = 0.0;
  if (!findRt(ctx, c, sqlite3_value_int(argv[0]), &rt_min)) {
    return;
  }
  double mass = sqlite3_value_double(argv[1]);
  int charge = sqlite3_value_int(argv[2]);
  if (charge <= 0) {
    sqlite3_result_error(ctx, "an envelope without a positive charge", -1);
    return;
  }
  sqlite3_result_double(
      ctx, c->model->correctMass(c->ms_level, mass, charge, rt_min));
}

// calib_mz_rt_ms(rt_ms, mz): the 3D tables, which keep the retention time in
// ms instead of a spectrum id
void calibMzRtMs(sqlite3_context* ctx, int, sqlite3_value** argv) {
  const LevelContext* c =
      static_cast<const LevelContext*>(sqlite3_user_data(ctx));
  double rt_min = sqlite3_value_double(argv[0]) / 60000.0;
  sqlite3_result_double(ctx, c->model->correctMz(
                                 c->ms_level, sqlite3_value_double(argv[1]),
                                 rt_min));
}

bool hasTable(sqlite3* db, const std::string& name) {
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db, "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?;");
  sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_TRANSIENT);
  bool found = sqlite3_step(stmt) == SQLITE_ROW;
  sqlite3_finalize(stmt);
  return found;
}

std::unordered_map<int, double> readRtMins(sqlite3* db,
                                           const std::string& table) {
  std::unordered_map<int, double> rt_mins;
  sqlite3_stmt* stmt =
      sql_util::prepareSql(db, "SELECT id, retention_time FROM " + table + ";");
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    rt_mins[sqlite3_column_int(stmt, 0)] = sqlite3_column_double(stmt, 1) / 60;
  }
  sqlite3_finalize(stmt);
  return rt_mins;
}

void registerFunctions(sqlite3* db, const std::string& suffix,
                       LevelContext* c) {
  int rc = sqlite3_create_function(db, ("calib_mz" + suffix).c_str(), 2,
                                   SQLITE_UTF8, c, calibMz, nullptr, nullptr);
  if (rc == SQLITE_OK) {
    rc = sqlite3_create_function(db, ("calib_mass" + suffix).c_str(), 3,
                                 SQLITE_UTF8, c, calibMass, nullptr, nullptr);
  }
  if (rc == SQLITE_OK) {
    rc = sqlite3_create_function(db, ("calib_mz_rt_ms" + suffix).c_str(), 2,
                                 SQLITE_UTF8, c, calibMzRtMs, nullptr,
                                 nullptr);
  }
  if (rc != SQLITE_OK) {
    LOG_ERROR("Cannot register the calibration functions: "
              << sqlite3_errmsg(db));
    std::exit(EXIT_FAILURE);
  }
}

// runs an UPDATE and returns the number of rows it changed
int update(sqlite3* db, const std::string& sql) {
  sql_util::execSql(db, sql);
  return sqlite3_changes(db);
}

void removeFile(const std::string& name) {
  std::error_code ec;
  std::filesystem::remove(name, ec);
}

}  // namespace

void process(const MassCalibModel& model, const std::string& ori_db_name,
             const std::string& calib_db_name) {
  auto start = std::chrono::steady_clock::now();
  std::string tmp_name = calib_db_name + ".tmp";
  for (const char* suffix : {"", "-journal", "-wal", "-shm"}) {
    removeFile(tmp_name + suffix);
  }
  // the copy goes through SQLite's backup API, not a file copy: TopFD
  // writes its database in WAL mode, and pages still in the -wal file are
  // part of it
  sqlite3* ori_db = nullptr;
  if (sqlite3_open_v2(ori_db_name.c_str(), &ori_db, SQLITE_OPEN_READONLY,
                      nullptr) != SQLITE_OK) {
    LOG_ERROR("Cannot open the database " << ori_db_name << ": "
                                          << sqlite3_errmsg(ori_db));
    std::exit(EXIT_FAILURE);
  }
  sqlite3* db = nullptr;
  if (sqlite3_open(tmp_name.c_str(), &db) != SQLITE_OK) {
    LOG_ERROR("Cannot open the database " << tmp_name << ": "
                                          << sqlite3_errmsg(db));
    std::exit(EXIT_FAILURE);
  }
  sqlite3_backup* backup = sqlite3_backup_init(db, "main", ori_db, "main");
  if (backup != nullptr) {
    sqlite3_backup_step(backup, -1);
    sqlite3_backup_finish(backup);
  }
  if (backup == nullptr || sqlite3_errcode(db) != SQLITE_OK) {
    LOG_ERROR("Cannot copy " << ori_db_name << " to " << tmp_name << ": "
                             << sqlite3_errmsg(db));
    std::exit(EXIT_FAILURE);
  }
  sqlite3_close(ori_db);
  // the file is a new copy under a temporary name: no journal is needed
  sql_util::execSql(db, "PRAGMA journal_mode = OFF;");
  sql_util::execSql(db, "PRAGMA synchronous = OFF;");

  // toppic's tables of an earlier run on TopFD's database; their pages stay
  // free in the file (no VACUUM: it may renumber the rows of the tables
  // without an INTEGER PRIMARY KEY, every TopFD table but the spectra)
  for (const std::string& table : kToppicTables) {
    if (hasTable(db, table)) {
      sql_util::execSql(db, "DROP TABLE " + table + ";");
    }
  }

  LevelContext ms1_ctx;
  ms1_ctx.model = &model;
  ms1_ctx.ms_level = 1;
  ms1_ctx.rt_mins = readRtMins(db, "ms1_spectrum");
  LevelContext ms2_ctx;
  ms2_ctx.model = &model;
  ms2_ctx.ms_level = 2;
  ms2_ctx.rt_mins = readRtMins(db, "ms2_spectrum");
  registerFunctions(db, "1", &ms1_ctx);
  registerFunctions(db, "2", &ms2_ctx);

  std::vector<std::pair<std::string, int>> counts;
  sql_util::execSql(db, "BEGIN;");
  for (int level = 1; level <= 2; level++) {
    std::string l = std::to_string(level);
    std::string p = "ms" + l + "_";
    counts.emplace_back(
        p + "peak", update(db, "UPDATE " + p + "peak SET mz = calib_mz" + l +
                                   "(spec_id, mz);"));
    counts.emplace_back(
        p + "env", update(db, "UPDATE " + p + "env SET mono_mass = calib_mass" +
                                  l + "(spec_id, mono_mass, charge), "
                                  "ref_mass = calib_mass" + l +
                                  "(spec_id, ref_mass, charge);"));
    counts.emplace_back(
        p + "env_peak", update(db, "UPDATE " + p + "env_peak SET mz = calib_mz" +
                                       l + "(spec_id, mz);"));
  }
  // the 3D tables of topfd -D: layer k is PEAKS<k>, its range the CONFIG row
  // inserted k-th (ms_sql_writer)
  bool has_config = hasTable(db, "CONFIG");
  for (int layer = 0; hasTable(db, "PEAKS" + std::to_string(layer));
       layer++) {
    std::string table = "PEAKS" + std::to_string(layer);
    counts.emplace_back(
        table, update(db, "UPDATE " + table +
                              " SET MZ = calib_mz_rt_ms1(RETENTIONTIME, MZ);"));
    if (has_config) {
      sql_util::execSql(db, "UPDATE CONFIG SET MZMIN = (SELECT MIN(MZ) FROM " +
                                table + "), MZMAX = (SELECT MAX(MZ) FROM " +
                                table + ") WHERE rowid = " +
                                std::to_string(layer + 1) + ";");
    }
  }

  sql_util::execSql(db,
                    "CREATE TABLE mass_calib(model TEXT NOT NULL,"
                    "source_database TEXT NOT NULL,"
                    "toppic_version TEXT NOT NULL,"
                    "time TEXT NOT NULL);");
  sqlite3_stmt* stmt = sql_util::prepareSql(
      db,
      "INSERT INTO mass_calib(model, source_database, toppic_version, time) "
      "VALUES (?, ?, ?, ?);");
  std::string json = model.toJson();
  std::string version = Version::getVersion();
  std::string time = time_util::getTimeStr();
  sqlite3_bind_text(stmt, 1, json.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 2, ori_db_name.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 3, version.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 4, time.c_str(), -1, SQLITE_TRANSIENT);
  sql_util::stepAndReset(db, stmt);
  sqlite3_finalize(stmt);
  sql_util::execSql(db, "COMMIT;");
  if (sqlite3_close(db) != SQLITE_OK) {
    LOG_ERROR("Cannot close the database " << tmp_name << ".");
    std::exit(EXIT_FAILURE);
  }
  std::error_code ec;
  std::filesystem::rename(tmp_name, calib_db_name, ec);
  if (ec) {
    LOG_ERROR("Cannot rename " << tmp_name << " to " << calib_db_name << ": "
                               << ec.message());
    std::exit(EXIT_FAILURE);
  }

  double seconds = std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - start)
                       .count();
  std::cout << "Mass calibration - database " << calib_db_name
            << " written from " << ori_db_name << " (";
  for (size_t i = 0; i < counts.size(); i++) {
    std::cout << (i == 0 ? "" : ", ") << counts[i].first << " "
              << counts[i].second;
  }
  std::cout << " rows corrected; " << std::fixed << std::setprecision(1)
            << seconds << " s)." << std::defaultfloat << std::endl;
}

}  // namespace sqlite_mass_calib

}  // namespace toppic
