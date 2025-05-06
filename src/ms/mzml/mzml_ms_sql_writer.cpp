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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/util/sql_util.hpp"
#include "ms/mzml/mzml_ms_sql_writer.hpp" 

namespace toppic {

namespace mzml_ms_sql_writer {

void writeMs2(sqlite3* sql_db, MzmlMsPtr ms_ptr, MatchEnvPtrVec &envs) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  //int ms_level = header_ptr->getMsLevel();
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

  std::string sql =
      "INSERT INTO ms2_spectrum(id, scan, retention_time, target_mz, begin_mz, end_mz, n_ion_type, c_ion_type) values ('" 
            + std::to_string(spec_id) + "',"
      + "'" + std::to_string(scan_num) + "',"
      + "'" + std::to_string(retention_time) + "',"
      + "'" + std::to_string(target_mz) + "',"
      + "'" + std::to_string(begin_mz) + "',"
      + "'" + std::to_string(end_mz) + "',"
      + "'" + n_ion_type + "',"
      + "'" + c_ion_type + "');";
  LOG_ERROR("INSERT SQL: " << sql); 
  sql_util::execSql(sql_db, sql); 

  PeakPtrVec raw_peaks = ms_ptr->getPeakPtrVec();
  for (size_t i = 0; i < raw_peaks.size(); i++) {
    std::string pos_str =
        str_util::fixedToString(raw_peaks[i]->getPosition(), 4);
    // peak.AddMember("mz", pos, allocator);
    std::string inte_str =
        str_util::toScientificStr(raw_peaks[i]->getIntensity(), 4);
    // peak.AddMember("intensity", inte, allocator);
    std::string sql =
        "INSERT INTO ms2_peak(spec_id, peak_id, mz, intensity) values ('" 
        + std::to_string(spec_id) + "'," 
        + "'" + std::to_string(i) + "'," 
        + "'" + pos_str + "'," 
        + "'" + inte_str + "');";
    sql_util::execSql(sql_db, sql);
  }

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

}

}
