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

#include "ms/mzml/mzml_ms_sql_writer.hpp"

#include <sqlite3.h>

#include <cstddef>
#include <cstdlib>
#include <mutex>
#include <string>

#include "common/util/logger.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace {

// Run a fully-bound INSERT and reset it so the handle can be reused.
void stepAndReset(sqlite3_stmt* stmt) {
  sqlite3_step(stmt);
  sqlite3_clear_bindings(stmt);
  sqlite3_reset(stmt);
}

}  // namespace

MzmlMsSqlWriter::MzmlMsSqlWriter(sqlite3* sql_db) : sql_db_(sql_db) {
  // Bulk-load PRAGMAs. WAL + synchronous=NORMAL removes the per-commit fsync
  // while staying crash-safe; the cache/mmap settings keep working pages in
  // memory. (This is a regenerable visualization database, so synchronous=OFF
  // with journal_mode=MEMORY would be faster still if durability is not
  // needed.)
  sql_util::execSql(sql_db_, "PRAGMA journal_mode = WAL;");
  sql_util::execSql(sql_db_, "PRAGMA synchronous = NORMAL;");
  sql_util::execSql(sql_db_, "PRAGMA temp_store = MEMORY;");
  sql_util::execSql(sql_db_, "PRAGMA cache_size = -65536;");    // ~64 MiB
  sql_util::execSql(sql_db_, "PRAGMA mmap_size = 268435456;");  // 256 MiB

  ms1_spec_stmt_ = sql_util::prepareSql(sql_db_,
      "INSERT INTO ms1_spectrum(id, scan, retention_time, peak_num, env_num, "
      "base_inte, min_ref_inte) VALUES (?, ?, ?, ?, ?, ?, ?);");
  ms1_peak_stmt_ = sql_util::prepareSql(sql_db_,
                           "INSERT INTO ms1_peak(spec_id, peak_id, mz, "
                           "intensity) VALUES (?, ?, ?, ?);");
  ms1_env_stmt_ = sql_util::prepareSql(sql_db_,
      "INSERT INTO ms1_env(spec_id, env_id, mono_mass, ref_mass, charge, "
      "intensity, envcnn_score, peak_num) VALUES (?, ?, ?, ?, ?, ?, ?, ?);");
  ms1_env_peak_stmt_ = sql_util::prepareSql(sql_db_,
      "INSERT INTO ms1_env_peak(spec_id, env_id, peak_id, mz, intensity) "
      "VALUES (?, ?, ?, ?, ?);");
  ms2_spec_stmt_ = sql_util::prepareSql(sql_db_,
      "INSERT INTO ms2_spectrum(id, scan, retention_time, target_mz, begin_mz, "
      "end_mz, n_ion_type, c_ion_type, peak_num, ms1_id) "
      "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);");
  ms2_peak_stmt_ = sql_util::prepareSql(sql_db_,
                           "INSERT INTO ms2_peak(spec_id, peak_id, mz, "
                           "intensity) VALUES (?, ?, ?, ?);");
  ms2_env_stmt_ = sql_util::prepareSql(sql_db_,
      "INSERT INTO ms2_env(spec_id, env_id, mono_mass, ref_mass, charge, "
      "intensity, envcnn_score, peak_num) VALUES (?, ?, ?, ?, ?, ?, ?, ?);");
  ms2_env_peak_stmt_ = sql_util::prepareSql(sql_db_,
      "INSERT INTO ms2_env_peak(spec_id, env_id, peak_id, mz, intensity) "
      "VALUES (?, ?, ?, ?, ?);");
}

MzmlMsSqlWriter::~MzmlMsSqlWriter() {
  flush();
  sqlite3_finalize(ms1_spec_stmt_);
  sqlite3_finalize(ms1_peak_stmt_);
  sqlite3_finalize(ms1_env_stmt_);
  sqlite3_finalize(ms1_env_peak_stmt_);
  sqlite3_finalize(ms2_spec_stmt_);
  sqlite3_finalize(ms2_peak_stmt_);
  sqlite3_finalize(ms2_env_stmt_);
  sqlite3_finalize(ms2_env_peak_stmt_);
}

void MzmlMsSqlWriter::begin() {
  if (!in_transaction_) {
    sql_util::execSql(sql_db_, "BEGIN TRANSACTION;");
    in_transaction_ = true;
  }
}

void MzmlMsSqlWriter::commit() {
  if (in_transaction_) {
    sql_util::execSql(sql_db_, "END TRANSACTION;");
    in_transaction_ = false;
    pending_ = 0;
  }
}

void MzmlMsSqlWriter::flush() {
  std::lock_guard<std::mutex> lock(mutex_);
  commit();
}

void MzmlMsSqlWriter::writeMs1(const MzmlMsPtr& ms_ptr,
                               const MatchEnvPtrVec& envs, double base_inte,
                               double min_ref_inte) {
  std::lock_guard<std::mutex> lock(mutex_);
  begin();

  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  int spec_id = header_ptr->getSpecId();
  const PeakPtrVec& raw_peaks = ms_ptr->getPeakPtrVec();

  sqlite3_bind_int(ms1_spec_stmt_, 1, spec_id);
  sqlite3_bind_int(ms1_spec_stmt_, 2, header_ptr->getFirstScanNum());
  sqlite3_bind_double(ms1_spec_stmt_, 3, header_ptr->getRetentionTime());
  sqlite3_bind_int(ms1_spec_stmt_, 4, static_cast<int>(raw_peaks.size()));
  sqlite3_bind_int(ms1_spec_stmt_, 5, static_cast<int>(envs.size()));
  sqlite3_bind_double(ms1_spec_stmt_, 6, base_inte);
  sqlite3_bind_double(ms1_spec_stmt_, 7, min_ref_inte);
  stepAndReset(ms1_spec_stmt_);

  for (size_t i = 0; i < raw_peaks.size(); i++) {
    sqlite3_bind_int(ms1_peak_stmt_, 1, spec_id);
    sqlite3_bind_int(ms1_peak_stmt_, 2, static_cast<int>(i));
    sqlite3_bind_double(ms1_peak_stmt_, 3, raw_peaks[i]->getPosition());
    sqlite3_bind_double(ms1_peak_stmt_, 4, raw_peaks[i]->getIntensity());
    stepAndReset(ms1_peak_stmt_);
  }

  for (size_t i = 0; i < envs.size(); i++) {
    EnvPtr theo_env = envs[i]->getTheoEnvPtr();
    int peak_num = theo_env->getPeakNum();
    sqlite3_bind_int(ms1_env_stmt_, 1, spec_id);
    sqlite3_bind_int(ms1_env_stmt_, 2, static_cast<int>(i));
    sqlite3_bind_double(ms1_env_stmt_, 3, theo_env->getMonoNeutralMass());
    // Neutral mass of the reference (most abundant) isotopic peak.
    sqlite3_bind_double(ms1_env_stmt_, 4, theo_env->getReferNeutralMass());
    sqlite3_bind_int(ms1_env_stmt_, 5, theo_env->getCharge());
    sqlite3_bind_double(ms1_env_stmt_, 6, theo_env->compInteSum());
    sqlite3_bind_double(ms1_env_stmt_, 7, envs[i]->getEnvcnnScore());
    sqlite3_bind_int(ms1_env_stmt_, 8, peak_num);
    stepAndReset(ms1_env_stmt_);

    for (int k = 0; k < peak_num; k++) {
      sqlite3_bind_int(ms1_env_peak_stmt_, 1, spec_id);
      sqlite3_bind_int(ms1_env_peak_stmt_, 2, static_cast<int>(i));
      sqlite3_bind_int(ms1_env_peak_stmt_, 3, k);
      sqlite3_bind_double(ms1_env_peak_stmt_, 4, theo_env->getMz(k));
      sqlite3_bind_double(ms1_env_peak_stmt_, 5, theo_env->getInte(k));
      stepAndReset(ms1_env_peak_stmt_);
    }
  }

  if (++pending_ >= COMMIT_CHUNK) {
    commit();
  }
}

void MzmlMsSqlWriter::writeMs2(const MzmlMsPtr& ms_ptr,
                               const MatchEnvPtrVec& envs) {
  std::lock_guard<std::mutex> lock(mutex_);
  begin();

  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  int spec_id = header_ptr->getSpecId();
  int ms1_id = header_ptr->getMsOneId();
  const PeakPtrVec& raw_peaks = ms_ptr->getPeakPtrVec();
  std::string n_ion_type =
      header_ptr->getActivationPtr()->getNIonTypePtr()->getName();
  std::string c_ion_type =
      header_ptr->getActivationPtr()->getCIonTypePtr()->getName();

  sqlite3_bind_int(ms2_spec_stmt_, 1, spec_id);
  sqlite3_bind_int(ms2_spec_stmt_, 2, header_ptr->getFirstScanNum());
  sqlite3_bind_double(ms2_spec_stmt_, 3, header_ptr->getRetentionTime());
  sqlite3_bind_double(ms2_spec_stmt_, 4, header_ptr->getPrecTargetMz());
  sqlite3_bind_double(ms2_spec_stmt_, 5, header_ptr->getPrecWinBegin());
  sqlite3_bind_double(ms2_spec_stmt_, 6, header_ptr->getPrecWinEnd());
  sqlite3_bind_text(ms2_spec_stmt_, 7, n_ion_type.c_str(), -1,
                    SQLITE_TRANSIENT);
  sqlite3_bind_text(ms2_spec_stmt_, 8, c_ion_type.c_str(), -1,
                    SQLITE_TRANSIENT);
  sqlite3_bind_int(ms2_spec_stmt_, 9, static_cast<int>(raw_peaks.size()));
  sqlite3_bind_int(ms2_spec_stmt_, 10, ms1_id);
  stepAndReset(ms2_spec_stmt_);

  for (size_t i = 0; i < raw_peaks.size(); i++) {
    sqlite3_bind_int(ms2_peak_stmt_, 1, spec_id);
    sqlite3_bind_int(ms2_peak_stmt_, 2, static_cast<int>(i));
    sqlite3_bind_double(ms2_peak_stmt_, 3, raw_peaks[i]->getPosition());
    sqlite3_bind_double(ms2_peak_stmt_, 4, raw_peaks[i]->getIntensity());
    stepAndReset(ms2_peak_stmt_);
  }

  for (size_t i = 0; i < envs.size(); i++) {
    EnvPtr theo_env = envs[i]->getTheoEnvPtr();
    int peak_num = theo_env->getPeakNum();
    sqlite3_bind_int(ms2_env_stmt_, 1, spec_id);
    sqlite3_bind_int(ms2_env_stmt_, 2, static_cast<int>(i));
    sqlite3_bind_double(ms2_env_stmt_, 3, theo_env->getMonoNeutralMass());
    // Neutral mass of the reference (most abundant) isotopic peak.
    sqlite3_bind_double(ms2_env_stmt_, 4, theo_env->getReferNeutralMass());
    sqlite3_bind_int(ms2_env_stmt_, 5, theo_env->getCharge());
    sqlite3_bind_double(ms2_env_stmt_, 6, theo_env->compInteSum());
    sqlite3_bind_double(ms2_env_stmt_, 7, envs[i]->getEnvcnnScore());
    sqlite3_bind_int(ms2_env_stmt_, 8, peak_num);
    stepAndReset(ms2_env_stmt_);

    for (int k = 0; k < peak_num; k++) {
      sqlite3_bind_int(ms2_env_peak_stmt_, 1, spec_id);
      sqlite3_bind_int(ms2_env_peak_stmt_, 2, static_cast<int>(i));
      sqlite3_bind_int(ms2_env_peak_stmt_, 3, k);
      sqlite3_bind_double(ms2_env_peak_stmt_, 4, theo_env->getMz(k));
      sqlite3_bind_double(ms2_env_peak_stmt_, 5, theo_env->getInte(k));
      stepAndReset(ms2_env_peak_stmt_);
    }
  }

  if (++pending_ >= COMMIT_CHUNK) {
    commit();
  }
}

}  // namespace toppic
