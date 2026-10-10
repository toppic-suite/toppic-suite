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

#ifndef TOPPIC_SQL_SQL_DB_COPY_HPP_
#define TOPPIC_SQL_SQL_DB_COPY_HPP_

#include <sqlite3.h>

#include <string>

namespace toppic {

namespace sql_db_copy {

// The copies of TopFD's database that toppic writes for its own spectrum
// files and searches in place of TopFD's (toppic_test experiment 154):
// <base>_mono.sqlite for the monoisotopic-mass copy of -g, and the
// corrected <base>_calib.sqlite (<base>_mono_calib.sqlite) of --mass-calib.
// TopFD's database is only read.

// Opens a copy of src_name, written under the temporary name
// dest_name + ".tmp" (an earlier temporary file is removed), for changes
// before closeCopy. The copy goes through SQLite's backup API, not a file
// copy: TopFD writes its database in WAL mode, and pages still in the -wal
// file are part of it. toppic's own tables of an earlier run (prsm,
// proteoform, protein, prsm_mass_shift, prsm_protein_match, fasta_seq;
// written again by every run) are dropped, without a VACUUM, which may
// renumber the rows of the tables without an INTEGER PRIMARY KEY. The copy
// has no journal. Exits with an error message when a database cannot be
// read or written.
sqlite3* openCopy(const std::string& src_name, const std::string& dest_name);

// Closes a copy opened by openCopy and renames it to dest_name, replacing an
// earlier file, so that an interrupted run leaves no partial file under that
// name.
void closeCopy(sqlite3* db, const std::string& dest_name);

// openCopy and closeCopy with no change in between.
void copy(const std::string& src_name, const std::string& dest_name);

}  // namespace sql_db_copy

}  // namespace toppic

#endif  // TOPPIC_SQL_SQL_DB_COPY_HPP_
