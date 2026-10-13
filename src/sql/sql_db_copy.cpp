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

#include "sql/sql_db_copy.hpp"

#include <cstdlib>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

#include "common/util/logger.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace sql_db_copy {

namespace {

// toppic's own tables, which every run writes again (prsm_sql_output)
const std::vector<std::string> kToppicTables = {
    "prsm",            "proteoform",         "protein",
    "prsm_mass_shift", "prsm_protein_match", "fasta_seq"};

std::string tmpName(const std::string& dest_name) {
  return dest_name + ".tmp";
}

void removeFile(const std::string& name) {
  std::error_code ec;
  std::filesystem::remove(name, ec);
}

}  // namespace

sqlite3* openCopy(const std::string& src_name, const std::string& dest_name) {
  std::string tmp_name = tmpName(dest_name);
  for (const char* suffix : {"", "-journal", "-wal", "-shm"}) {
    removeFile(tmp_name + suffix);
  }
  // read only as immutable: no -wal / -shm files next to TopFD's database
  sqlite3* src_db = nullptr;
  if (sqlite3_open_v2(sql_util::immutableUri(src_name).c_str(), &src_db,
                      SQLITE_OPEN_READONLY | SQLITE_OPEN_URI,
                      nullptr) != SQLITE_OK) {
    LOG_ERROR("Cannot open the database " << src_name << ": "
                                          << sqlite3_errmsg(src_db));
    std::exit(EXIT_FAILURE);
  }
  sqlite3* db = nullptr;
  if (sqlite3_open(tmp_name.c_str(), &db) != SQLITE_OK) {
    LOG_ERROR("Cannot open the database " << tmp_name << ": "
                                          << sqlite3_errmsg(db));
    std::exit(EXIT_FAILURE);
  }
  sqlite3_backup* backup = sqlite3_backup_init(db, "main", src_db, "main");
  if (backup != nullptr) {
    sqlite3_backup_step(backup, -1);
    sqlite3_backup_finish(backup);
  }
  if (backup == nullptr || sqlite3_errcode(db) != SQLITE_OK) {
    LOG_ERROR("Cannot copy " << src_name << " to " << tmp_name << ": "
                             << sqlite3_errmsg(db));
    std::exit(EXIT_FAILURE);
  }
  sqlite3_close(src_db);
  // the file is a new copy under a temporary name: no journal is needed
  sql_util::execSql(db, "PRAGMA journal_mode = OFF;");
  sql_util::execSql(db, "PRAGMA synchronous = OFF;");
  // their pages stay free in the file
  for (const std::string& table : kToppicTables) {
    if (sql_util::hasTable(db, table)) {
      sql_util::execSql(db, "DROP TABLE " + table + ";");
    }
  }
  return db;
}

void closeCopy(sqlite3* db, const std::string& dest_name) {
  std::string tmp_name = tmpName(dest_name);
  // WAL mode, as TopFD writes its database: a reader of the copy (a viewer)
  // then does not block toppic's later writes to it (toppic_test experiment
  // 156: in rollback mode the post mass matching stopped with "database is
  // locked")
  sql_util::execSql(db, "PRAGMA journal_mode = WAL;");
  if (sqlite3_close(db) != SQLITE_OK) {
    LOG_ERROR("Cannot close the database " << tmp_name << ".");
    std::exit(EXIT_FAILURE);
  }
  std::error_code ec;
  std::filesystem::rename(tmp_name, dest_name, ec);
  if (ec) {
    LOG_ERROR("Cannot rename " << tmp_name << " to " << dest_name << ": "
                               << ec.message());
    std::exit(EXIT_FAILURE);
  }
}

void copyTables(const std::string& src_name, const std::string& dest_name,
                const std::vector<std::string>& tables) {
  std::string tmp_name = tmpName(dest_name);
  for (const char* suffix : {"", "-journal", "-wal", "-shm"}) {
    removeFile(tmp_name + suffix);
  }
  sqlite3* db = nullptr;
  if (sqlite3_open_v2(
          tmp_name.c_str(), &db,
          SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_URI,
          nullptr) != SQLITE_OK) {
    LOG_ERROR("Cannot open the database " << tmp_name << ": "
                                          << sqlite3_errmsg(db));
    std::exit(EXIT_FAILURE);
  }
  // a new file under a temporary name: no journal is needed
  sql_util::execSql(db, "PRAGMA journal_mode = OFF;");
  sql_util::execSql(db, "PRAGMA synchronous = OFF;");
  sqlite3_stmt* stmt = sql_util::prepareSql(db, "ATTACH DATABASE ? AS src;");
  std::string src_uri = sql_util::immutableUri(src_name);
  sqlite3_bind_text(stmt, 1, src_uri.c_str(), -1, SQLITE_TRANSIENT);
  if (sqlite3_step(stmt) != SQLITE_DONE) {
    LOG_ERROR("Cannot open the database " << src_name << ": "
                                          << sqlite3_errmsg(db));
    std::exit(EXIT_FAILURE);
  }
  sqlite3_finalize(stmt);
  sql_util::execSql(db, "BEGIN;");
  for (const std::string& table : tables) {
    stmt = sql_util::prepareSql(
        db,
        "SELECT sql FROM src.sqlite_master WHERE type = 'table' AND name = ?;");
    sqlite3_bind_text(stmt, 1, table.c_str(), -1, SQLITE_TRANSIENT);
    if (sqlite3_step(stmt) != SQLITE_ROW) {
      LOG_ERROR("The database " << src_name << " has no table " << table
                                << ".");
      std::exit(EXIT_FAILURE);
    }
    std::string create_sql =
        reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
    sqlite3_finalize(stmt);
    // the statement names no schema: the table is created in main
    sql_util::execSql(db, create_sql + ";");
    sql_util::execSql(db, "INSERT INTO main." + table + " SELECT * FROM src." +
                              table + " ORDER BY rowid;");
  }
  sql_util::execSql(db, "COMMIT;");
  sql_util::execSql(db, "DETACH DATABASE src;");
  closeCopy(db, dest_name);
}

bool removeDb(const std::string& name) {
  bool found = std::filesystem::exists(name);
  for (const char* suffix : {"", "-journal", "-wal", "-shm"}) {
    removeFile(name + suffix);
  }
  return found;
}

}  // namespace sql_db_copy

}  // namespace toppic
