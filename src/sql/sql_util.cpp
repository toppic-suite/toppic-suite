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

#include "sql/sql_util.hpp"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

#include "common/util/logger.hpp"

namespace toppic {

namespace sql_util {

void execSql(sqlite3* sql_db, const std::string& sql) {
  char* errMsg = 0;
  int rc;
  // Execute SQL statement
  rc = sqlite3_exec(sql_db, sql.c_str(), 0, 0, &errMsg);
  if (rc != SQLITE_OK) {
    LOG_ERROR("Sql" << sql);
    LOG_ERROR("SQL error: " << errMsg);
    sqlite3_free(errMsg);
    exit(EXIT_FAILURE);
  }
}

sqlite3_stmt* prepareSql(sqlite3* sql_db, const std::string& sql) {
  sqlite3_stmt* stmt = nullptr;
  if (sqlite3_prepare_v2(sql_db, sql.c_str(), -1, &stmt, nullptr) !=
      SQLITE_OK) {
    LOG_ERROR("Failed to prepare SQL: " << sql);
    LOG_ERROR("SQL error: " << sqlite3_errmsg(sql_db));
    exit(EXIT_FAILURE);
  }
  return stmt;
}

void stepAndReset(sqlite3* sql_db, sqlite3_stmt* stmt) {
  if (sqlite3_step(stmt) != SQLITE_DONE) {
    LOG_ERROR("SQL error: " << sqlite3_errmsg(sql_db));
    exit(EXIT_FAILURE);
  }
  sqlite3_reset(stmt);
}

bool hasTable(sqlite3* sql_db, const std::string& name) {
  sqlite3_stmt* stmt = prepareSql(
      sql_db, "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?;");
  sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_TRANSIENT);
  bool found = sqlite3_step(stmt) == SQLITE_ROW;
  sqlite3_finalize(stmt);
  return found;
}

std::string immutableUri(const std::string& file_name) {
  std::string path = std::filesystem::path(file_name).generic_string();
  std::string uri;
  for (char c : path) {
    if (c == '%' || c == '?' || c == '#') {
      char buf[4];
      std::snprintf(buf, sizeof(buf), "%%%02X", static_cast<unsigned char>(c));
      uri += buf;
    } else {
      uri += c;
    }
  }
  // a Windows drive letter needs an empty authority: file:///C:/...
  if (uri.size() >= 2 && std::isalpha(static_cast<unsigned char>(uri[0])) &&
      uri[1] == ':') {
    uri = "///" + uri;
  }
  return "file:" + uri + "?immutable=1";
}

}  // namespace sql_util

}  // namespace toppic
