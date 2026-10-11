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

#ifndef TOPPIC_PRSM_PRSM_SQL_OUTPUT_HPP_
#define TOPPIC_PRSM_PRSM_SQL_OUTPUT_HPP_

#include <string>

#include "para/prsm_para.hpp"

namespace toppic {

namespace prsm_sql_output {

// Writes the identifications of a search into the SQLite database
// sql_file_name (toppic's <base>_toppic.sqlite of the spectrum file,
// prsm_post_mass_match::toppicSqlFileName; topmg's: TopFD's database,
// prsm_post_mass_match::sqlFileName): the prsm table (with mass-shift and
// protein-match details) from <base_name>.<prsm_file_ext>, the proteoform
// table from <base_name>.<form_file_ext>, the protein table from
// <base_name>.<prot_file_ext>, and the fasta_seq table from the FASTA file
// (base_name is the spectrum file name without its extension). Earlier
// tables of the same names are replaced. Does nothing (with a message) when
// sql_file_name is empty (a combined result of -c) or the database does not
// exist, e.g. when topfd was run with --no-sql.
void write(const PrsmParaPtr& prsm_para_ptr, const std::string& sql_file_name,
           const std::string& base_name, const std::string& fasta_file_name,
           const std::string& prsm_file_ext,
           const std::string& form_file_ext,
           const std::string& prot_file_ext);

}  // namespace prsm_sql_output

}  // namespace toppic

#endif
