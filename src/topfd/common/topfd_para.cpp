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

#include <iostream>
#include <sstream>
#include <iomanip>
#include "common/util/version.hpp"
#include "common/util/time_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/util/sql_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/version.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

TopfdPara::~TopfdPara() {
  if (sql_db != nullptr) {
    sqlite3_close(sql_db);
  }
}

void TopfdPara::setMzmlFileNameAndFaims(std::string &mzml_file_name, 
                                        bool is_faims, double voltage) {
  mzml_file_name_ = mzml_file_name;
  is_faims_ = is_faims;
  faims_volt_ = voltage;
  output_base_name_ = file_util::basename(mzml_file_name_);
  // if it is faims data, then add integer voltage to output_file_name
  if (is_faims_) {
    output_base_name_ = output_base_name_ + "_" 
      + str_util::toString(static_cast<int>(faims_volt_)); 
  }
  html_dir_ =  output_base_name_ + "_" + "html";
  ms1_json_dir_ = html_dir_ 
    + file_util::getFileSeparator() + "topfd" 
    + file_util::getFileSeparator() + "ms1_json";
  ms2_json_dir_ = html_dir_ 
    + file_util::getFileSeparator() + "topfd" 
    + file_util::getFileSeparator() + "ms2_json";
  sql_file_name_ = output_base_name_ + ".sqlite";
  createSqlDb(sql_file_name_);
}

void TopfdPara::createSqlDb(std::string sql_db_name) {
  int rc;
  // Open database
  if (sql_db != nullptr) {
    sqlite3_close(sql_db);
  }
  rc = sqlite3_open(sql_db_name.c_str(), &sql_db);
  if (rc) {
    LOG_ERROR("Can't open database: " << sqlite3_errmsg(sql_db)); 
    exit(EXIT_FAILURE);
  }

  std::string sql = "CREATE TABLE IF NOT EXISTS ms1_spectrum(id INTEGER PRIMARY KEY,"
                                                            "scan INTEGER NOT NULL,"
                                                            "retention_time REAL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db, sql); 
  sql = "CREATE TABLE IF NOT EXISTS ms1_peak(spec_id INTEGER NOT NULL,"
                                            "peak_id INTEGER NOT NULL,"
                                            "mz REAL NOT NULL,"
                                            "intensity REAL NOT NULL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db, sql);
  sql = "DELETE from ms1_peak;"; 
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db, sql);
  sql = "DELETE from ms1_spectrum;"; 
  LOG_DEBUG("SQL: " << sql);

  sql = "CREATE TABLE IF NOT EXISTS ms2_spectrum(id INTEGER PRIMARY KEY,"
        "scan INTEGER NOT NULL,"
        "retention_time REAL,"
        "target_mz REAL,"
        "begin_mz REAL,"
        "end_mz REAL,"
        "n_ion_type TEXT,"
        "c_ion_type TEXT);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db, sql); 
  sql = "CREATE TABLE IF NOT EXISTS ms2_peak(spec_id INTEGER NOT NULL,"
                                            "peak_id INTEGER NOT NULL,"
                                            "mz REAL NOT NULL,"
                                            "intensity REAL NOT NULL);";
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db, sql);
  sql = "DELETE from ms2_peak;"; 
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db, sql);
  sql = "DELETE from ms2_spectrum;"; 
  LOG_DEBUG("SQL: " << sql);
  sql_util::execSql(sql_db, sql); 
}

std::string TopfdPara::getTopfdParaStr(const std::string &prefix,
                                       const std::string &sep,
                                       int gap) {
  std::stringstream output;
  output << prefix << std::setw(gap) << std::left 
      << "File name:                  " << sep  << mzml_file_name_ << std::endl;
  if (is_faims_) {
    output << prefix << std::setw(gap) << std::left 
      << "Faims data:                 " << sep << "Yes" << std::endl;
    output << prefix << std::setw(gap) << std::left 
      << "Faims voltage:              " << sep << faims_volt_ << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Faims data:                 " << sep << "No" << std::endl;
    output << prefix << std::setw(gap) << std::left 
      << "Faims voltage:              " << sep << "N/A"<< std::endl;
  }
  output << prefix << std::setw(gap) << std::left 
      << "Number of MS1 scans:        " << sep  << ms_1_scan_num_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Number of MS/MS scans:      " << sep  << ms_2_scan_num_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Spectral data type:         " << sep  << "Centroid" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Maximum charge:             "  << sep << max_charge_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Maximum monoisotopic mass:  " << sep << max_mass_ << " Dalton" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Peak m/z error tolerance:   " << sep << mz_error_ << " m/z" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Thread number:              " << sep << thread_num_ << std::endl;
  if (gene_html_folder_) {
    output << prefix << std::setw(gap) << std::left 
      << "Generate Html files:        " << sep << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Generate Html files:        " << sep << "No" << std::endl;
  }

  if (missing_level_one_) {
    output << prefix << std::setw(gap) << std::left 
      << "Miss MS1 spectra:           " << sep << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Miss MS1 spectra:           " << sep << "No" << std::endl;

    output << std::endl << prefix << std::setw(gap) << std::left
      << "####### MS1 spectral deconvolution parameters   #######" << std::endl;
    output << prefix << std::setw(gap) << std::left 
      << "MS1 signal/noise ratio:     " << sep << ms_one_sn_ratio_ << std::endl;
    output << prefix << std::setw(gap) << std::left
      << "####### MS1 spectral deconvolution parameters   #######" << std::endl;

    output << std::endl << prefix << std::setw(gap) << std::left
      << "####### MS1 feature detection parameters        #######" << std::endl;
    output << prefix << std::setw(gap) << std::left 
      << "Feature min scan number:                        " << sep << ms1_min_scan_num_ << std::endl;

    if (use_single_scan_noise_level_) {
      output << prefix << std::setw(gap) << std::left 
      << "Use single scan noise level:                    " << sep << "Yes" << std::endl;
    }
    else {
      output << prefix << std::setw(gap) << std::left 
      << "Use single scan noise level:                    " << sep << "No" << std::endl;
    }
    output << prefix << std::setw(gap) << std::left
      << "Intensity ratio for splitting features:         " << sep << split_intensity_ratio_ << std::endl;

    output << prefix << std::setw(gap) << std::left 
      << "Feature ECScore cutoff:                         " << sep  << ms1_ecscore_cutoff_ << std::endl;
    if (search_prec_window_) {
      output << prefix << std::setw(gap) << std::left 
      << "Additional feature search for isolation windows:" << sep << "Yes" << std::endl;
    }
    else {
      output << prefix << std::setw(gap) << std::left 
      << "Additional feature search for isolation windows:" << sep << "No" << std::endl;
    }
    output << prefix << std::setw(gap) << std::left
      << "####### MS1 feature detection parameters        #######" << std::endl;
  }

  output << std::endl << prefix << std::setw(gap) << std::left
      << "####### MS/MS spectral deconvolution parameters #######" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Default precursor window:                           " << sep << prec_window_ << " m/z" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Activation type:                                    " << sep  << activation_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "MS/MS signal/noise ratio:                           " << sep << ms_two_sn_ratio_ << std::endl;
  if (sort_use_msdeconv_) {
    output << prefix << std::setw(gap) << std::left 
      << "Fragment envelope ranking:                          " << sep << "MS-Deconv score" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Fragment envelope ranking:                          " << sep << "EnvCNN score" << std::endl;
  }
  output << prefix << std::setw(gap) << std::left 
      << "Fragment envelope EnvCNN score cutoff:              " << sep  << ms2_env_cnn_score_cutoff_<< std::endl;
  if (aa_num_based_filter_) {
    output << prefix << std::setw(gap) << std::left 
      << "Filtering fragments using estimated fragment number:" << sep << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Filtering fragments using estimated fragment number:" << sep << "No" << std::endl;
  }
  output << prefix << std::setw(gap) << std::left
      << "####### MS/MS spectral deconvolution parameters #######" << std::endl;

  return output.str();
}

std::string TopfdPara::getParaStr(const std::string &prefix, 
		                  const std::string &sep) {
  std::stringstream output;
  int gap = 25;
  output << prefix << "TopFD " << Version::getVersion() << std::endl;
  output << prefix << "Timestamp: " << time_util::getTimeStr() << std::endl;
  output << prefix << "###################### Parameters #####################" << std::endl;
  output << getTopfdParaStr(prefix, sep, gap);
  output << prefix << "###################### Parameters #####################" << std::endl;
  return output.str();
}

}  // namespace toppic
