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

#include "topfd/common/topfd_para.hpp"

#include <sqlite3.h>

#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/time_util.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

// A banner line with the title centered and padded with '#' to a fixed width,
// e.g. "############### Parameters ###############".
std::string TopfdPara::banner(const std::string& prefix,
                              const std::string& title) {
  int fill = para_banner_width_ - 2 - static_cast<int>(title.size());
  if (fill < 2) fill = 2;
  int left = fill / 2;
  int right = fill - left;
  return prefix + std::string(left, '#') + " " + title + " " +
         std::string(right, '#');
}

TopfdPara::~TopfdPara() {
  if (sql_db_ != nullptr) {
    sqlite3_close(sql_db_);
  }
}

void TopfdPara::setMzmlFileNameAndFaims(const std::string& mzml_file_name,
                                        bool is_faims, double voltage) {
  mzml_file_name_ = mzml_file_name;
  is_faims_ = is_faims;
  faims_volt_ = voltage;
  output_base_name_ = file_util::basename(mzml_file_name_);
  // if it is faims data, then add integer voltage to output_file_name
  if (is_faims_) {
    output_base_name_ =
        output_base_name_ + "_" + std::to_string(static_cast<int>(faims_volt_));
  }
  sql_file_name_ = output_base_name_ + ".sqlite";
  if (gene_sql_) {
    createSqlDb(sql_file_name_);
  }
}

void TopfdPara::createSqlDb(const std::string& sql_db_name) {
  int rc;
  // Open database
  if (sql_db_ != nullptr) {
    sqlite3_close(sql_db_);
  }
  if (std::filesystem::exists(sql_db_name)) {
    file_util::delFile(sql_db_name);
  }
  rc = sqlite3_open(sql_db_name.c_str(), &sql_db_);
  if (rc) {
    LOG_ERROR("Can't open database: " << sqlite3_errmsg(sql_db_));
    exit(EXIT_FAILURE);
  }

  std::string sql =
      "CREATE TABLE IF NOT EXISTS ms1_spectrum(id INTEGER PRIMARY KEY,"
      "scan INTEGER NOT NULL,"
      "retention_time REAL,"
      "peak_num INTEGER,"
      "env_num INTEGER,"
      "base_inte REAL,"
      "min_ref_inte REAL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);
  sql =
      "CREATE TABLE IF NOT EXISTS ms1_peak(spec_id INTEGER NOT NULL,"
      "peak_id INTEGER NOT NULL,"
      "mz REAL NOT NULL,"
      "intensity REAL NOT NULL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);
  sql = "DELETE from ms1_peak;";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);
  sql = "DELETE from ms1_spectrum;";
  LOG_DEBUG("SQL: " << sql);

  sql =
      "CREATE TABLE IF NOT EXISTS ms2_spectrum(id INTEGER PRIMARY KEY,"
      "scan INTEGER NOT NULL,"
      "retention_time REAL,"
      "target_mz REAL,"
      "begin_mz REAL,"
      "end_mz REAL,"
      "n_ion_type TEXT,"
      "c_ion_type TEXT,"
      "peak_num INTEGER,"
      "ms1_id INTEGER,"
      "FOREIGN KEY(ms1_id) REFERENCES ms1_spectrum(id));";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);
  sql =
      "CREATE TABLE IF NOT EXISTS ms2_peak(spec_id INTEGER NOT NULL,"
      "peak_id INTEGER NOT NULL,"
      "mz REAL NOT NULL,"
      "intensity REAL NOT NULL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);
  sql =
      "CREATE TABLE IF NOT EXISTS ms_info(ms1_scan_num INTEGER NOT NULL,"
      "ms2_scan_num INTEGER NOT NULL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);
  sql = "INSERT INTO ms_info(ms1_scan_num, ms2_scan_num) values ('" +
        std::to_string(ms_1_scan_num_) + "'," + "'" +
        std::to_string(ms_2_scan_num_) + "');";
  LOG_DEBUG("INSERT SQL: " << sql);
  sql_util::execSql(sql_db_, sql);

  sql =
      "CREATE TABLE IF NOT EXISTS ms1_env(spec_id INTEGER NOT NULL,"
      "env_id INTEGER NOT NULL,"
      "mono_mass REAL NOT NULL,"
      "charge INTEGER NOT NULL,"
      "intensity REAL NOT NULL,"
      "envcnn_score REAL NOT NULL,"
      "peak_num INTEGER NOT NULL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);

  sql =
      "CREATE TABLE IF NOT EXISTS ms1_env_peak(spec_id INTEGER NOT NULL,"
      "env_id INTEGER NOT NULL,"
      "peak_id INTEGER NOT NULL,"
      "mz REAL NOT NULL,"
      "intensity REAL NOT NULL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);

  sql =
      "CREATE TABLE IF NOT EXISTS ms2_env(spec_id INTEGER NOT NULL,"
      "env_id INTEGER NOT NULL,"
      "mono_mass REAL NOT NULL,"
      "charge INTEGER NOT NULL,"
      "intensity REAL NOT NULL,"
      "envcnn_score REAL NOT NULL,"
      "peak_num INTEGER NOT NULL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);

  sql =
      "CREATE TABLE IF NOT EXISTS ms2_env_peak(spec_id INTEGER NOT NULL,"
      "env_id INTEGER NOT NULL,"
      "peak_id INTEGER NOT NULL,"
      "mz REAL NOT NULL,"
      "intensity REAL NOT NULL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db_, sql);
}

std::string TopfdPara::getTopfdParaStr(const std::string& prefix,
                                       const std::string& sep) const {
  std::stringstream output;
  const int w = para_label_width_;
  auto kv = [&](const char* label) -> std::ostream& {
    return output << prefix << std::setw(w) << std::left << label << sep;
  };

  kv("File name:") << mzml_file_name_ << std::endl;
  if (is_faims_) {
    kv("FAIMS data:") << "Yes" << std::endl;
    kv("FAIMS voltage:") << faims_volt_ << std::endl;
  } else {
    kv("FAIMS data:") << "No" << std::endl;
    kv("FAIMS voltage:") << "N/A" << std::endl;
  }
  kv("Number of MS1 scans:") << ms_1_scan_num_ << std::endl;
  kv("Number of MS/MS scans:") << ms_2_scan_num_ << std::endl;
  kv("Spectral data type:") << "Centroid" << std::endl;
  kv("Maximum charge:") << max_charge_ << std::endl;
  kv("Maximum monoisotopic mass:") << max_mass_ << " Dalton" << std::endl;
  kv("Peak m/z error tolerance:") << mz_error_ << " m/z" << std::endl;
  kv("Thread number:") << thread_num_ << std::endl;

  if (missing_level_one_) {
    kv("Miss MS1 spectra:") << "Yes" << std::endl;
  } else {
    kv("Miss MS1 spectra:") << "No" << std::endl;

    output << std::endl
           << banner(prefix, "MS1 spectral deconvolution parameters")
           << std::endl;
    kv("MS1 signal/noise ratio:") << ms_one_sn_ratio_ << std::endl;

    output << std::endl
           << banner(prefix, "MS1 feature detection parameters") << std::endl;
    kv("Feature min scan number:") << ms1_min_scan_num_ << std::endl;
    kv("Use single scan noise level:")
        << (use_single_scan_noise_level_ ? "Yes" : "No") << std::endl;
    kv("Intensity ratio for splitting features:")
        << split_intensity_ratio_ << std::endl;
    kv("Feature ECScore cutoff:") << ms1_ecscore_cutoff_ << std::endl;
    kv("Additional feature search for isolation windows:")
        << (search_prec_window_ ? "Yes" : "No") << std::endl;
  }

  output << std::endl
         << banner(prefix, "MS/MS spectral deconvolution parameters")
         << std::endl;
  if (isFilePrecWindow()) {
    // the input file carries the MS/MS precursor windows, so prec_window_ (the
    // default width) is not used.
    kv("Default precursor window:") << "FILE" << std::endl;
  } else {
    kv("Default precursor window:") << prec_window_ << " m/z" << std::endl;
  }
  kv("Activation type:") << activation_ << std::endl;
  kv("MS/MS signal/noise ratio:") << ms_two_sn_ratio_ << std::endl;
  kv("Fragment envelope ranking:")
      << (sort_use_msdeconv_ ? "MS-Deconv score" : "EnvCNN score") << std::endl;
  kv("Fragment envelope EnvCNN score cutoff:")
      << ms2_env_cnn_score_cutoff_ << std::endl;
  kv("Filtering fragments using estimated fragment number:")
      << (aa_num_based_filter_ ? "Yes" : "No") << std::endl;

  return output.str();
}

std::string TopfdPara::getParaStr(const std::string& prefix,
                                  const std::string& sep) const {
  std::stringstream output;
  output << prefix << "Timestamp: " << time_util::getTimeStr() << std::endl;
  output << banner(prefix, "Parameters") << std::endl;
  output << getTopfdParaStr(prefix, sep);
  output << banner(prefix, "Parameters") << std::endl;
  return output.str();
}

}  // namespace toppic
