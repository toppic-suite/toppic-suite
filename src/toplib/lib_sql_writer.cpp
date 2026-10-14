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

#include "toplib/lib_builder.hpp"

#include <cstdlib>
#include <filesystem>
#include <string>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

#include <sqlite3.h>

#include "common/util/logger.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace toplib {

namespace lib_sql_writer {

namespace {

const char kTableName[] = "target_decoy_spectra_representatives";

const char* typeName(LibColType t) {
  switch (t) {
    case LibColType::kInteger:
      return "INTEGER";
    case LibColType::kReal:
      return "REAL";
    default:
      return "TEXT";
  }
}

// pandas quotes identifiers with double quotes
std::string quote(const std::string& name) { return "\"" + name + "\""; }

void bindValue(sqlite3_stmt* stmt, int pos, const LibValue& v, LibColType t) {
  if (std::holds_alternative<std::monostate>(v)) {
    sqlite3_bind_null(stmt, pos);
  } else if (const long long* l = std::get_if<long long>(&v)) {
    if (t == LibColType::kInteger) {
      sqlite3_bind_int64(stmt, pos, *l);
    } else {
      sqlite3_bind_double(stmt, pos, static_cast<double>(*l));
    }
  } else if (const double* d = std::get_if<double>(&v)) {
    sqlite3_bind_double(stmt, pos, *d);
  } else {
    const std::string& s = std::get<std::string>(v);
    sqlite3_bind_text(stmt, pos, s.c_str(), static_cast<int>(s.size()),
                      SQLITE_TRANSIENT);
  }
}

}  // namespace

void write(const LibLibrary& lib, const std::string& db_file_name) {
  // the columns: representative_id, the precursor, 50 masses, intensities,
  // charges and normalized intensities, the target/decoy flag, lib_ptm
  std::vector<std::pair<std::string, LibColType>> cols = {
      {"representative_id", LibColType::kInteger},
      {"precursor_mass", LibColType::kReal},
      {"precursor_intensity", LibColType::kReal},
      {"precursor_charge", LibColType::kInteger}};
  for (int k = 0; k < kLibPeakNum; k++) {
    cols.push_back({"mass_" + std::to_string(k), LibColType::kReal});
  }
  for (int k = 0; k < kLibPeakNum; k++) {
    cols.push_back({"intensity_" + std::to_string(k), LibColType::kReal});
  }
  // a charge column with a missing value is float64 in pandas
  for (int k = 0; k < kLibPeakNum; k++) {
    bool full = true;
    for (const LibRep& rep : lib.reps) {
      full &= rep.charges.size() > static_cast<size_t>(k);
    }
    cols.push_back({"charge_" + std::to_string(k),
                    full ? LibColType::kInteger : LibColType::kReal});
  }
  for (int k = 0; k < kLibPeakNum; k++) {
    cols.push_back({"norm_intensity_" + std::to_string(k), LibColType::kReal});
  }
  cols.push_back({"flag", LibColType::kInteger});
  for (const LibIdColumn& c : lib.ptm_columns) {
    cols.push_back({c.name, c.type});
  }

  std::filesystem::path parent = std::filesystem::path(db_file_name).parent_path();
  if (!parent.empty()) {
    std::error_code ec;
    std::filesystem::create_directories(parent, ec);
    if (ec) {
      LOG_ERROR("Cannot create the directory " << parent.string() << ": "
                                               << ec.message());
      exit(EXIT_FAILURE);
    }
  }
  sqlite3* db = nullptr;
  if (sqlite3_open(db_file_name.c_str(), &db) != SQLITE_OK) {
    LOG_ERROR("Cannot open the database " << db_file_name << ": "
                                          << sqlite3_errmsg(db));
    sqlite3_close(db);
    exit(EXIT_FAILURE);
  }
  // pandas' to_sql(if_exists='replace'): drop the table, then create it
  sql_util::execSql(db, "DROP TABLE IF EXISTS " + quote(kTableName));
  std::string create = "CREATE TABLE " + quote(kTableName) + " (\n";
  std::string insert = "INSERT INTO " + quote(kTableName) + " (";
  std::string params;
  for (size_t i = 0; i < cols.size(); i++) {
    if (i > 0) {
      create += ",\n  ";
      insert += ",";
      params += ",";
    }
    create += quote(cols[i].first) + " " + typeName(cols[i].second);
    insert += quote(cols[i].first);
    params += "?";
  }
  create += "\n)";
  insert += ") VALUES (" + params + ")";
  sql_util::execSql(db, create);

  sql_util::execSql(db, "BEGIN TRANSACTION");
  sqlite3_stmt* stmt = sql_util::prepareSql(db, insert);
  // the targets, then the decoys, of the identified representatives
  for (int flag = 1; flag >= 0; flag--) {
    for (size_t r = 0; r < lib.reps.size(); r++) {
      const LibRep& rep = lib.reps[r];
      if (!rep.identified) {
        continue;
      }
      int pos = 1;
      sqlite3_bind_int64(stmt, pos++, static_cast<long long>(r) + 1);
      sqlite3_bind_double(stmt, pos++, rep.precursor_mass);
      sqlite3_bind_double(stmt, pos++, rep.precursor_intensity);
      sqlite3_bind_int64(stmt, pos++, rep.precursor_charge);
      const std::vector<double>& masses = flag == 1 ? rep.masses : rep.decoy_masses;
      size_t n = rep.masses.size();
      for (size_t k = 0; k < static_cast<size_t>(kLibPeakNum); k++, pos++) {
        if (k < n) {
          sqlite3_bind_double(stmt, pos, masses[k]);
        } else {
          sqlite3_bind_null(stmt, pos);
        }
      }
      for (size_t k = 0; k < static_cast<size_t>(kLibPeakNum); k++, pos++) {
        if (k < n) {
          sqlite3_bind_double(stmt, pos, rep.intensities[k]);
        } else {
          sqlite3_bind_null(stmt, pos);
        }
      }
      for (size_t k = 0; k < static_cast<size_t>(kLibPeakNum); k++, pos++) {
        if (k >= n) {
          sqlite3_bind_null(stmt, pos);
        } else if (cols[pos - 1].second == LibColType::kInteger) {
          sqlite3_bind_int64(stmt, pos, rep.charges[k]);
        } else {
          sqlite3_bind_double(stmt, pos, rep.charges[k]);
        }
      }
      for (size_t k = 0; k < static_cast<size_t>(kLibPeakNum); k++, pos++) {
        if (k < n) {
          sqlite3_bind_double(stmt, pos, rep.norm_intensities[k]);
        } else {
          sqlite3_bind_null(stmt, pos);
        }
      }
      sqlite3_bind_int64(stmt, pos++, flag);
      for (size_t c = 0; c < rep.ptm_values.size(); c++, pos++) {
        bindValue(stmt, pos, rep.ptm_values[c], lib.ptm_columns[c].type);
      }
      sql_util::stepAndReset(db, stmt);
    }
  }
  sqlite3_finalize(stmt);
  sql_util::execSql(db, "COMMIT");
  sqlite3_close(db);
}

}  // namespace lib_sql_writer

}  // namespace toplib

}  // namespace toppic
