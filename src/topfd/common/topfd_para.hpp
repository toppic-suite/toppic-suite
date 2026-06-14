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

#ifndef TOPPIC_TOPFD_COMMON_TOPFD_PARA_HPP_
#define TOPPIC_TOPFD_COMMON_TOPFD_PARA_HPP_

#include <sqlite3.h>

#include <memory>
#include <string>
#include <utility>

namespace toppic {

class TopfdPara;

using TopfdParaPtr = std::shared_ptr<TopfdPara>;

class TopfdPara {
 public:
  TopfdPara() = default;

  ~TopfdPara();

  std::string getTopfdParaStr(const std::string& prefix,
                              const std::string& sep) const;

  std::string getParaStr(const std::string& prefix,
                         const std::string& sep) const;

  void createSqlDb(const std::string& db_name);

  const std::string& getExeDir() const { return exe_dir_; }
  const std::string& getResourceDir() const { return resource_dir_; }
  bool isMissingLevelOne() const { return missing_level_one_; }
  int getMaxCharge() const { return max_charge_; }
  double getMaxMass() const { return max_mass_; }
  double getMzError() const { return mz_error_; }

  bool isEstimateMinInte() const { return estimate_min_inte_; }
  double getMsOneSnRatio() const { return ms_one_sn_ratio_; }
  double getMsTwoSnRatio() const { return ms_two_sn_ratio_; }
  double getSplitIntensityRatio() const { return split_intensity_ratio_; }
  double getPrecWindowWidth() const { return prec_window_; }
  bool isSortUseMsDeconv() const { return sort_use_msdeconv_; }
  bool isAANumBasedFilter() const { return aa_num_based_filter_; }
  double getMs2EnvCnnScoreCutoff() const { return ms2_env_cnn_score_cutoff_; }
  const std::string& getActivation() const { return activation_; }
  bool isGeneSql() const { return gene_sql_; }
  bool isKeepUnusedPeaks() const { return keep_unused_peaks_; }
  bool isOutputMultipleMass() const { return output_multiple_mass_; }
  bool isOutputCsvFeatureFile() const { return output_csv_feature_file_; }
  bool isOutputDpEnvs() const { return output_dp_envs_; }
  int getThreadNum() const { return thread_num_; }
  double getMs1EcscoreCutoff() const { return ms1_ecscore_cutoff_; }
  double getMs2EcscoreCutoff() const { return ms2_ecscore_cutoff_; }
  bool isSearchPrecWindow() const { return search_prec_window_; }
  // True when the input mzML carries MS/MS precursor windows (so prec_window_,
  // the default width, is unused). Set from MzmlProfile::hasPrecWindow().
  bool isFilePrecWindow() const { return file_prec_window_; }
  bool isUseSingleScanNoiseLevel() const {
    return use_single_scan_noise_level_;
  }
  bool isTextPeakList() const { return text_peak_list_; }
  double getPrecInteCutoffRatio() const { return prec_inte_cutoff_ratio_; }

  const std::string& getMzmlFileName() const { return mzml_file_name_; }
  const std::string& getOutputBaseName() const { return output_base_name_; }

  int getFracId() const { return frac_id_; }
  bool isFaims() const { return is_faims_; }
  double getFaimsVoltage() const { return faims_volt_; }
  int getMs1ScanNum() const { return ms_1_scan_num_; }
  int getMs2ScanNum() const { return ms_2_scan_num_; }
  int getMs1MinScanNum() const { return ms1_min_scan_num_; }
  int getMs2MinScanNum() const { return ms2_min_scan_num_; }
  sqlite3* getSqlDb() const { return sql_db_; }

  void setExeDir(std::string dir) { exe_dir_ = std::move(dir); }
  void setResourceDir(std::string dir) { resource_dir_ = std::move(dir); }
  void setMissingLevelOne(bool missing) { missing_level_one_ = missing; }
  void setMaxCharge(int charge) { max_charge_ = charge; }
  void setMaxMass(double mass) { max_mass_ = mass; }
  void setMzError(double error) { mz_error_ = error; }
  void setMsOneSnRatio(double ratio) { ms_one_sn_ratio_ = ratio; }
  void setMsTwoSnRatio(double ratio) { ms_two_sn_ratio_ = ratio; }
  void setSplitIntensityRatio(double ratio) { split_intensity_ratio_ = ratio; }
  void setPrecWindowWidth(double window) { prec_window_ = window; }
  void setSortUseMsDeconv(bool use) { sort_use_msdeconv_ = use; }
  void setAANumBasedFilter(bool filter) { aa_num_based_filter_ = filter; }
  void setMs2EnvCnnScoreCutoff(double cutoff) {
    ms2_env_cnn_score_cutoff_ = cutoff;
  }
  void setActivation(std::string activation) {
    activation_ = std::move(activation);
  }
  void setKeepUnusedPeaks(bool keep) { keep_unused_peaks_ = keep; }
  void setOutputMultipleMass(bool output) { output_multiple_mass_ = output; }
  void setOutputCsvFeatureFile(bool output) {
    output_csv_feature_file_ = output;
  }
  void setGeneSql(bool gene_sql) { gene_sql_ = gene_sql; }
  void setOutputDpEnvs(bool output) { output_dp_envs_ = output; }
  void setThreadNum(int num) { thread_num_ = num; }
  void setSearchPrecWindow(bool search) { search_prec_window_ = search; }
  void setFilePrecWindow(bool file_prec_window) {
    file_prec_window_ = file_prec_window;
  }
  void setUseSingleScanNoiseLevel(bool single_scan_noise) {
    use_single_scan_noise_level_ = single_scan_noise;
  }
  void setMs1EcscoreCutoff(double cutoff) { ms1_ecscore_cutoff_ = cutoff; }
  void setMs2EcscoreCutoff(double cutoff) { ms2_ecscore_cutoff_ = cutoff; }
  void setMs1MinScanNum(int min_scan_num) { ms1_min_scan_num_ = min_scan_num; }
  void setMs2MinScanNum(int min_scan_num) { ms2_min_scan_num_ = min_scan_num; }

  void setFracId(int frac_id) { frac_id_ = frac_id; }
  void setMzmlFileNameAndFaims(const std::string& mzml_file_name, bool is_faims,
                               double voltage);
  void setMs1ScanNumber(int ms1_scan_num) { ms_1_scan_num_ = ms1_scan_num; }
  void setMs2ScanNumber(int ms2_scan_num) { ms_2_scan_num_ = ms2_scan_num; }

  void setTextPeakList(bool text_peak_list) {
    text_peak_list_ = text_peak_list;
  }
  void setOutputMatchEnv(bool output_match_env) {
    output_match_env_ = output_match_env;
  }

 private:
  // Width of the label column in the parameter printout: every value starts at
  // the same column so the report is aligned across all sections. It is wide
  // enough for the longest label ("Filtering fragments using estimated fragment
  // number:").
  static constexpr int para_label_width_ = 53;

  // Total width of the "### <title> ###" section banners.
  static constexpr int para_banner_width_ = 55;

  // A banner line with the title centered and padded with '#' to a fixed width,
  // e.g. "############### Parameters ###############".
  static std::string banner(const std::string& prefix,
                            const std::string& title);

  std::string exe_dir_;
  std::string resource_dir_;

  // parameters for deconcovolution
  int max_charge_ = 30;
  double max_mass_ = 50000;
  // precursor window is used only when the mzML file does not
  // contain the precursor window information
  double prec_window_ = 3.0;
  bool missing_level_one_ = false;
  double mz_error_ = 0.02;
  double ms_one_sn_ratio_ = 3.0;
  double ms_two_sn_ratio_ = 1.0;
  int thread_num_ = 1;
  std::string activation_ = "FILE";
  // keep unused peaks in dynamic programming
  bool keep_unused_peaks_ = false;
  // sorting using msdeconv, the default method is env_cnn score
  bool sort_use_msdeconv_ = false;
  double ms2_env_cnn_score_cutoff_ = 0.0;
  bool aa_num_based_filter_ = true;
  bool output_csv_feature_file_ = false;
  bool gene_sql_ = false;
  // set per input file from MzmlProfile::hasPrecWindow()
  bool file_prec_window_ = false;

  // parameters for feature identification
  double split_intensity_ratio_ = 2.5;
  bool search_prec_window_ = true;
  bool use_single_scan_noise_level_ = false;
  double ms1_ecscore_cutoff_ = 0.1;
  double ms2_ecscore_cutoff_ = 0;
  int ms1_min_scan_num_ = 1;
  int ms2_min_scan_num_ = 1;

  // For an MS/MS spectrum, the precursor is not reported if its intensity
  // is less than the cutoff ratio * the intensity of the first precusor
  double prec_inte_cutoff_ratio_ = 0.1;

  //** Fixed parameter setting **
  // estimate min intensity using the method in Thrash.
  bool estimate_min_inte_ = true;
  bool output_multiple_mass_ = false;
  bool output_match_env_ = false;
  // When true, dump the windowed candidate envelopes and the DP-selected
  // envelopes of every deconvoluted spectrum to win_envs.txt / dp_envs.txt
  // (debugging aid).
  bool output_dp_envs_ = false;

  //** information for each run **
  int frac_id_ = -1;
  std::string mzml_file_name_ = "";
  bool is_faims_ = false;
  double faims_volt_ = -1;
  std::string output_base_name_ = "deconv";
  std::string sql_file_name_ = "";
  sqlite3* sql_db_ = nullptr;

  int ms_1_scan_num_ = -1;
  int ms_2_scan_num_ = -1;

  // call function for processing a simple text file
  // containing a mass list
  bool text_peak_list_ = false;
};

}  // namespace toppic

#endif
