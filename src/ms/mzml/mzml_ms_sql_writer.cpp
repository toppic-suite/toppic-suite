//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include <boost/thread/mutex.hpp>
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "sql/sql_util.hpp"
#include "ms/mzml/mzml_ms_sql_writer.hpp" 

namespace toppic {

namespace mzml_ms_sql_writer {

// serialization mutex.
boost::mutex writer_mutex;

void writeMs1(sqlite3* sql_db, MzmlMsPtr ms_ptr, MatchEnvPtrVec& envs, double base_inte, double min_ref_inte) {
  char * err_msg = 0; 
  const char *tail_peak;
  sqlite3_stmt *stmt_peak;
  std::string sql_peak = "INSERT INTO ms1_peak(spec_id, peak_id, mz, intensity) VALUES (?, ?, ?, ?);";

  const char *tail_env;
  sqlite3_stmt *stmt_env;
  std::string sql_env = "INSERT INTO ms1_env(spec_id, env_id, mono_mass, charge, intensity, envcnn_score, peak_num) VALUES (?, ?, ?, ?, ?, ?, ?);";

  const char *tail_env_peak;
  sqlite3_stmt *stmt_env_peak;
  std::string sql_env_peak = "INSERT INTO ms1_env_peak(spec_id, env_id, peak_id, mz, intensity) VALUES (?, ?, ?, ?, ?);";

  writer_mutex.lock();
  sqlite3_prepare_v2(sql_db, sql_peak.c_str(), 256, &stmt_peak, &tail_peak);
  sqlite3_prepare_v2(sql_db, sql_env.c_str(), 256, &stmt_env, &tail_env);
  sqlite3_prepare_v2(sql_db, sql_env_peak.c_str(), 256, &stmt_env_peak, &tail_env_peak);

  sqlite3_exec(sql_db, "BEGIN TRANSACTION", NULL, NULL, &err_msg);

  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  int spec_id = header_ptr->getSpecId();
  int scan_num = header_ptr->getFirstScanNum();
  double retention_time = header_ptr->getRetentionTime();
  PeakPtrVec raw_peaks = ms_ptr->getPeakPtrVec();
  std::string sql = "INSERT INTO ms1_spectrum(id, scan, retention_time, peak_num, env_num, base_inte, min_ref_inte) values ('" 
  + std::to_string(spec_id) + "'," 
  + "'" + std::to_string(scan_num) + "'," 
  + "'" + std::to_string(retention_time) + "',"
  + "'" + std::to_string(raw_peaks.size()) + "',"
  + "'" + std::to_string(envs.size()) + "',"
  + "'" + std::to_string(base_inte) + "',"
  + "'" + std::to_string(min_ref_inte) + "');";
  LOG_DEBUG("INSERT SQL: " << sql);
  sql_util::execSql(sql_db, sql);

  for (size_t i = 0; i < raw_peaks.size(); i++) {
    sqlite3_bind_int(stmt_peak, 1, spec_id); 
    sqlite3_bind_int(stmt_peak, 2, i); 
    sqlite3_bind_double(stmt_peak, 3, raw_peaks[i]->getPosition()); 
    sqlite3_bind_double(stmt_peak, 4, raw_peaks[i]->getIntensity()); 
    sqlite3_step(stmt_peak);
    sqlite3_clear_bindings(stmt_peak);
    sqlite3_reset(stmt_peak);
  }

  for (size_t i = 0; i < envs.size(); i++) {
    EnvPtr theo_env = envs[i]->getTheoEnvPtr();
    sqlite3_bind_int(stmt_env, 1, spec_id); 
    sqlite3_bind_int(stmt_env, 2, i); 
    sqlite3_bind_double(stmt_env, 3, theo_env->getMonoNeutralMass());  
    sqlite3_bind_int(stmt_env, 4, theo_env->getCharge()); 
    sqlite3_bind_double(stmt_env, 5, theo_env->compInteSum());  
    sqlite3_bind_double(stmt_env, 6, envs[i]->getEnvcnnScore());  
    sqlite3_bind_double(stmt_env, 7, theo_env->getPeakNum());  
    sqlite3_step(stmt_env);
    sqlite3_clear_bindings(stmt_env);
    sqlite3_reset(stmt_env);
    for (int k = 0; k < theo_env->getPeakNum(); k++) {
      sqlite3_bind_int(stmt_env_peak, 1, spec_id);
      sqlite3_bind_int(stmt_env_peak, 2, i);
      sqlite3_bind_int(stmt_env_peak, 3, k); 
      sqlite3_bind_double(stmt_env_peak, 4, theo_env->getMz(k));
      sqlite3_bind_double(stmt_env_peak, 5, theo_env->getInte(k)); 
      sqlite3_step(stmt_env_peak);
      sqlite3_clear_bindings(stmt_env_peak);
      sqlite3_reset(stmt_env_peak);
    }
  }

  sqlite3_exec(sql_db, "END TRANSACTION", NULL, NULL, &err_msg);
  writer_mutex.unlock();
}

void writeMs2(sqlite3* sql_db, MzmlMsPtr ms_ptr, MatchEnvPtrVec &envs) {
  char * err_msg = 0; 
  const char *tail;
  sqlite3_stmt *stmt;
  std::string sql = "INSERT INTO ms2_peak(spec_id, peak_id, mz, intensity) VALUES (?, ?, ?, ?);";
  writer_mutex.lock();
  sqlite3_prepare_v2(sql_db, sql.c_str(), 256, &stmt, &tail);

  sqlite3_exec(sql_db, "BEGIN TRANSACTION", NULL, NULL, &err_msg);
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  int spec_id = header_ptr->getSpecId();
  int scan_num = header_ptr->getFirstScanNum();
  double retention_time = header_ptr->getRetentionTime();
  double target_mz = header_ptr->getPrecTargetMz();
  double begin_mz = header_ptr->getPrecWinBegin();
  double end_mz = header_ptr->getPrecWinEnd();
  std::string n_ion_type =
      header_ptr->getActivationPtr()->getNIonTypePtr()->getName();
  std::string c_ion_type =
      header_ptr->getActivationPtr()->getCIonTypePtr()->getName();
  PeakPtrVec raw_peaks = ms_ptr->getPeakPtrVec();

  sql =
      "INSERT INTO ms2_spectrum(id, scan, retention_time, target_mz, begin_mz, end_mz, n_ion_type, c_ion_type, peak_num) values ('" 
            + std::to_string(spec_id) + "',"
      + "'" + std::to_string(scan_num) + "',"
      + "'" + std::to_string(retention_time) + "',"
      + "'" + std::to_string(target_mz) + "',"
      + "'" + std::to_string(begin_mz) + "',"
      + "'" + std::to_string(end_mz) + "',"
      + "'" + n_ion_type + "',"
      + "'" + c_ion_type + "',"
      + "'" + std::to_string(raw_peaks.size()) +"');";
  LOG_DEBUG("INSERT SQL: " << sql); 
  sql_util::execSql(sql_db, sql); 

  for (size_t i = 0; i < raw_peaks.size(); i++) {
    sqlite3_bind_int(stmt, 1, spec_id); 
    sqlite3_bind_int(stmt, 2, i); 
    sqlite3_bind_double(stmt, 3, raw_peaks[i]->getPosition()); 
    sqlite3_bind_double(stmt, 4, raw_peaks[i]->getIntensity()); 
    sqlite3_step(stmt);
    sqlite3_clear_bindings(stmt);
    sqlite3_reset(stmt);
  }


  sqlite3_exec(sql_db, "END TRANSACTION", NULL, NULL, &err_msg);

  writer_mutex.unlock();

  /*
  rapidjson::Value envelopes(rapidjson::kArrayType);
  for (size_t i = 0; i < envs.size(); i++) {
    rapidjson::Value env(rapidjson::kObjectType);
    EnvPtr theo_env = envs[i]->getTheoEnvPtr();
    env.AddMember("id", i, allocator);
    env.AddMember("mono_mass", theo_env->getMonoNeutralMass(), allocator);
    env.AddMember("charge", theo_env->getCharge(), allocator);

    rapidjson::Value env_peaks(rapidjson::kArrayType);
    for (int k = 0; k < theo_env->getPeakNum(); k++) {
      rapidjson::Value peak(rapidjson::kObjectType);
      peak.AddMember("mz", theo_env->getMz(k), allocator);
      peak.AddMember("intensity", theo_env->getInte(k), allocator);
      env_peaks.PushBack(peak, allocator);
    }
    env.AddMember("env_peaks", env_peaks, allocator);
    envelopes.PushBack(env, allocator);
  }
  doc.AddMember("envelopes", envelopes, allocator);
  */
}

}  // namespace mzml_ms_sql_writer
}
