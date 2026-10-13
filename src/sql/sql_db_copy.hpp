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
#include <vector>

namespace toppic {

namespace sql_db_copy {

// toppic's and topmg's database of a spectrum file, <base>_toppic.sqlite and
// <base>_topmg.sqlite, built from TopFD's (toppic_test experiments 154, 156,
// 168 and 170): the tables ms2_env and ms2_env_peak (writeEnvelopeDb), to
// which toppic's post mass matching appends its envelopes, or with toppic
// --mass-calib the corrected copy of every table (openCopy and closeCopy,
// sqlite_mass_calib). TopFD's database is only read.

// Opens a copy of src_name, written under the temporary name
// dest_name + ".tmp" (an earlier temporary file is removed), for changes
// before closeCopy. The copy goes through SQLite's backup API, not a file
// copy; src_name is opened read only as immutable (no -wal / -shm files next
// to it; it must not change meanwhile). toppic's own tables of an earlier run
// (prsm, proteoform, protein, prsm_mass_shift, prsm_protein_match, fasta_seq;
// written again by every run) are dropped, without a VACUUM, which may
// renumber the rows of the tables without an INTEGER PRIMARY KEY. A
// compressed src_name (topfd's default, toppic_test experiment 172) is
// copied decompressed (zstd_db::decompressDb), with its peaks and envelopes
// as rows. The copy has no journal. Exits with an error message when a
// database cannot be read or written.
sqlite3* openCopy(const std::string& src_name, const std::string& dest_name);

// Closes a copy opened by openCopy in WAL mode (TopFD's journal mode, so that
// a reader of the copy does not block a later writer) and renames it to
// dest_name, replacing an earlier file, so that an interrupted run leaves no
// partial file under that name.
void closeCopy(sqlite3* db, const std::string& dest_name);

// Writes dest_name, a new database with the tables of src_name (rowid
// tables) and their rows in rowid order: each table is created with
// src_name's own CREATE statement, src_name is attached read only as
// immutable (no -wal / -shm files next to it; it must not change meanwhile),
// and the new database is written under the temporary name and renamed as
// closeCopy does, replacing an earlier file. A table that a compressed
// src_name (topfd's default, toppic_test experiment 172) stores as BLOB rows
// is written as rows again, with its own CREATE statement
// (zstd_db::BlobTableReader). Exits with an error message when a database
// cannot be read or written or src_name lacks a table.
void copyTables(const std::string& src_name, const std::string& dest_name,
                const std::vector<std::string>& tables);

// Writes dest_name with the tables ms2_env and ms2_env_peak of TopFD's
// database topfd_db_name (copyTables); when TopFD's database does not exist
// (topfd -z), prints a note and writes nothing. Returns whether dest_name was
// written.
bool writeEnvelopeDb(const std::string& topfd_db_name,
                     const std::string& dest_name);

// Deletes the database name with its -journal, -wal and -shm files, if they
// exist. Returns whether name existed.
bool removeDb(const std::string& name);

}  // namespace sql_db_copy

}  // namespace toppic

#endif  // TOPPIC_SQL_SQL_DB_COPY_HPP_
