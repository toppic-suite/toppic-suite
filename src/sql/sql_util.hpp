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

#ifndef TOPPIC_SQL_SQL_UTIL_HPP_
#define TOPPIC_SQL_SQL_UTIL_HPP_

#include <sqlite3.h>

#include <string>

namespace toppic {

namespace sql_util {

// Execute a statement, logging and aborting on failure.
void execSql(sqlite3* sql_db, const std::string& sql);

// Compile a statement, aborting on failure (matching execSql).
sqlite3_stmt* prepareSql(sqlite3* sql_db, const std::string& sql);

// Run a fully-bound statement (aborting on failure) and reset it so the
// handle can be reused; every parameter must be rebound before the next step.
void stepAndReset(sqlite3* sql_db, sqlite3_stmt* stmt);

}  // namespace sql_util

}  // namespace toppic

#endif
