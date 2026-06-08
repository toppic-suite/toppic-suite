//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_MS_MZML_MZML_MS_SQL_WRITER_HPP_
#define TOPPIC_MS_MZML_MZML_MS_SQL_WRITER_HPP_

#include <mutex>

#include <sqlite3.h>

#include "ms/env/match_env.hpp"
#include "ms/mzml/mzml_ms.hpp"

namespace toppic {

// Writes deconvoluted MS1/MS2 spectra (peaks and matched envelopes) into an
// SQLite database.
//
// Performance: rather than a transaction-per-spectrum (one fsync per scan), the
// spectra are inserted across a small number of large transactions using
// prepared statements that are compiled once and reused. Call writeMs1/writeMs2
// once per spectrum; the writer commits automatically every kCommitChunk
// spectra and flushes the final batch when it is destroyed (or via flush()).
//
// The sqlite3 connection is owned by the caller, which creates the schema before
// constructing the writer and closes the database after destroying it. The
// writer only manages the bulk-load PRAGMAs, prepared statements and the
// chunked transactions on that connection.
class MzmlMsSqlWriter {
 public:
  explicit MzmlMsSqlWriter(sqlite3 *sql_db);

  ~MzmlMsSqlWriter();

  MzmlMsSqlWriter(const MzmlMsSqlWriter &) = delete;
  MzmlMsSqlWriter &operator=(const MzmlMsSqlWriter &) = delete;

  void writeMs1(const MzmlMsPtr &ms_ptr, const MatchEnvPtrVec &envs,
                double base_inte, double min_ref_inte);

  void writeMs2(const MzmlMsPtr &ms_ptr, const MatchEnvPtrVec &envs);

  // Commit any pending spectra (also called by the destructor).
  void flush();

 private:
  void begin();   // open a transaction if none is open (caller holds mutex_)
  void commit();  // commit the open transaction, if any (caller holds mutex_)

  sqlite3 *sql_db_;   // not owned
  std::mutex mutex_;  // serializes writers sharing this object

  sqlite3_stmt *ms1_spec_stmt_ = nullptr;
  sqlite3_stmt *ms1_peak_stmt_ = nullptr;
  sqlite3_stmt *ms1_env_stmt_ = nullptr;
  sqlite3_stmt *ms1_env_peak_stmt_ = nullptr;
  sqlite3_stmt *ms2_spec_stmt_ = nullptr;
  sqlite3_stmt *ms2_peak_stmt_ = nullptr;

  int pending_ = 0;             // spectra inserted in the current transaction
  bool in_transaction_ = false;

  // Commit (and reopen) a transaction every this many spectra, to amortize the
  // fsync across many scans while bounding the WAL/journal size.
  static constexpr int kCommitChunk = 2000;
};

// Backward-compatible free-function facade matching the historic
// writeMs1/writeMs2(sqlite3*, ...) signature, for callers that pass a raw
// connection rather than holding a MzmlMsSqlWriter. One writer is kept per
// sqlite3 connection behind the scenes, so prepared statements and chunked
// transactions are still reused across calls (these are NOT per-call writers).
//
// IMPORTANT: because writes are batched, the caller MUST call close(sql_db)
// before sqlite3_close(sql_db). close() commits the final partial transaction
// and finalizes the prepared statements (sqlite3_close fails while statements
// are live). Skipping it loses the last < kCommitChunk spectra. Prefer the
// MzmlMsSqlWriter object directly in new code; it flushes in its destructor.
namespace mzml_ms_sql_writer {

void writeMs1(sqlite3 *sql_db, const MzmlMsPtr &ms_ptr, const MatchEnvPtrVec &envs,
              double base_inte, double min_ref_inte);

void writeMs2(sqlite3 *sql_db, const MzmlMsPtr &ms_ptr, const MatchEnvPtrVec &envs);

// Flush and release the writer for this connection. Call before sqlite3_close.
void close(sqlite3 *sql_db);

}  // namespace mzml_ms_sql_writer

}  // namespace toppic

#endif
