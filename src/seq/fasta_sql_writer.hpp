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

#ifndef TOPPIC_SEQ_FASTA_SQL_WRITER_HPP_
#define TOPPIC_SEQ_FASTA_SQL_WRITER_HPP_

#include <sqlite3.h>

#include <string>

namespace toppic {

namespace fasta_sql_writer {

// (Re)creates the fasta_seq table (name, description, sequence) in the
// database and fills it with the proteins of the FASTA file.
void write(sqlite3* sql_db, const std::string& fasta_file_name);

}  // namespace fasta_sql_writer

}  // namespace toppic

#endif
