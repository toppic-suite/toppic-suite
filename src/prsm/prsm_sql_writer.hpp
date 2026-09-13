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

#ifndef TOPPIC_PRSM_PRSM_SQL_WRITER_HPP_
#define TOPPIC_PRSM_PRSM_SQL_WRITER_HPP_

#include <sqlite3.h>

#include <memory>
#include <string>

#include "para/prsm_para.hpp"
#include "prsm/prsm.hpp"
#include "prsm/search_fasta_match.hpp"

namespace toppic {

// Writes identified PrSMs into an SQLite database (the topfd .sqlite file of
// the spectrum file), mirroring the columns of the PrSM/proteoform TSV tables
// written by PrsmMatchTableWriter. The sqlite3 connection is owned by the
// caller.
class PrsmSqlWriter {
 public:
  PrsmSqlWriter(const PrsmParaPtr& prsm_para_ptr, sqlite3* sql_db);

  // Reads the PrSMs of <spectrum base name>.<input_file_ext> and (re)creates
  // the table table_name with one row per PrSM. With write_details, the
  // prsm_mass_shift table (one row per mass shift of a PrSM) and the
  // prsm_protein_match table (the other proteins whose sequence also matches
  // the PrSM) are (re)created as well; the rows of both are keyed by prsm_id.
  void write(const std::string& input_file_ext, const std::string& table_name,
             bool write_details);

 private:
  void createPrsmTable(const std::string& table_name);
  void writePrsm(sqlite3_stmt* stmt, const PrsmPtr& prsm_ptr);
  void writeMassShifts(sqlite3_stmt* stmt, const PrsmPtr& prsm_ptr);
  void writeProteinMatches(sqlite3_stmt* stmt, const PrsmPtr& prsm_ptr);

  PrsmParaPtr prsm_para_ptr_;
  sqlite3* sql_db_;  // not owned
  SearchFastaMatchPtr search_match_ptr_;
};

using PrsmSqlWriterPtr = std::shared_ptr<PrsmSqlWriter>;

}  // namespace toppic

#endif
