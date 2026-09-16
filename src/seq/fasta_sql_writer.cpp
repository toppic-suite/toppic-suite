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

#include "seq/fasta_sql_writer.hpp"

#include <sqlite3.h>

#include <string>

#include "seq/fasta_reader.hpp"
#include "seq/fasta_seq.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace fasta_sql_writer {

void write(sqlite3* sql_db, const std::string& fasta_file_name) {
  sql_util::execSql(sql_db, "DROP TABLE IF EXISTS fasta_seq;");
  sql_util::execSql(sql_db,
                    "CREATE TABLE fasta_seq("
                    "name TEXT PRIMARY KEY,"
                    "description TEXT NOT NULL,"
                    "sequence TEXT NOT NULL);");
  sqlite3_stmt* stmt = sql_util::prepareSql(
      sql_db,
      "INSERT OR REPLACE INTO fasta_seq(name, description, sequence) "
      "VALUES (?, ?, ?);");
  sql_util::execSql(sql_db, "BEGIN;");
  FastaReader reader(fasta_file_name);
  for (FastaSeqPtr seq_ptr = reader.getNextSeq(); seq_ptr != nullptr;
       seq_ptr = reader.getNextSeq()) {
    sqlite3_bind_text(stmt, 1, seq_ptr->getName().c_str(), -1,
                      SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, seq_ptr->getDesc().c_str(), -1,
                      SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 3, seq_ptr->getRawSeq().c_str(), -1,
                      SQLITE_TRANSIENT);
    sql_util::stepAndReset(sql_db, stmt);
  }
  reader.close();
  sql_util::execSql(sql_db, "COMMIT;");
  sqlite3_finalize(stmt);
}

}  // namespace fasta_sql_writer

}  // namespace toppic
