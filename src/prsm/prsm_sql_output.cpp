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

#include "prsm/prsm_sql_output.hpp"

#include <sqlite3.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "para/prsm_para.hpp"
#include "prsm/prsm_sql_writer.hpp"
#include "seq/fasta_sql_writer.hpp"

namespace toppic {

namespace prsm_sql_output {

void write(const PrsmParaPtr& prsm_para_ptr, const std::string& sp_file_name,
           const std::string& fasta_file_name,
           const std::string& prsm_file_ext,
           const std::string& form_file_ext,
           const std::string& prot_file_ext) {
  std::string sql_base = file_util::basename(sp_file_name);
  // the spectra written by the post mass matching keep TopFD's database name
  if (str_util::endsWith(sql_base, "_post_ms2")) {
    sql_base = sql_base.substr(0, sql_base.size() - 9);
  } else if (str_util::endsWith(sql_base, "_ms2")) {
    sql_base = sql_base.substr(0, sql_base.size() - 4);
  }
  std::string sql_file_name = sql_base + ".sqlite";
  if (!std::filesystem::exists(sql_file_name)) {
    std::cout << "SQLite database " << sql_file_name
              << " not found: identifications are not written to a database."
              << std::endl;
    return;
  }
  std::cout << "Writing identifications to " << sql_file_name << " - started."
            << std::endl;
  sqlite3* sql_db = nullptr;
  if (sqlite3_open(sql_file_name.c_str(), &sql_db) != SQLITE_OK) {
    LOG_ERROR("Cannot open the database " << sql_file_name << ": "
                                          << sqlite3_errmsg(sql_db));
    exit(EXIT_FAILURE);
  }
  {
    PrsmSqlWriter sql_writer(prsm_para_ptr, sql_db);
    sql_writer.write(prsm_file_ext, "prsm", true);
    sql_writer.write(form_file_ext, "proteoform", false);
    sql_writer.write(prot_file_ext, "protein", false);
  }
  fasta_sql_writer::write(sql_db, fasta_file_name);
  sqlite3_close(sql_db);
  std::cout << "Writing identifications to " << sql_file_name << " - finished."
            << std::endl;
}

}  // namespace prsm_sql_output

}  // namespace toppic
