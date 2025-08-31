//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include "common/util/logger.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace sql_util {

void execSql(sqlite3 *sql_db, const std::string &sql) {
  char *errMsg = 0;
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

}  // namespace sql_util

}  // namespace toppic

