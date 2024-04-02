//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_TOPDIA_COMMON_TOPDIA_PARA_HPP_
#define TOPPIC_TOPDIA_COMMON_TOPDIA_PARA_HPP_

#include <memory>
#include <string>

namespace toppic {
class TopdiaPara {
public:
  TopdiaPara() {};

  std::string getParaStr(const std::string &prefix,
                         const std::string &sep);

  std::string getExeDir() {return exe_dir_;}
  std::string getResourceDir() {return resource_dir_;}
  bool isMissingLevelOne() {return missing_level_one_;}
  int getMaxCharge() {return max_charge_;}
  double getMaxMass() {return max_mass_;}
  double getMzError() {return mz_error_;}

  bool isEstimateMinInte() {return estimate_min_inte_;}
  double getMsOneSnRatio() {return ms_one_sn_ratio_;}
  double getMsTwoSnRatio() {return ms_two_sn_ratio_;}
  double getPrecWindowWidth() {return prec_window_;}
  bool isUseMsDeconv() {return use_msdeconv_;}
  bool isDoFinalFiltering() {return do_final_filtering_;}
  std::string getActivation() {return activation_;}
  bool isGeneHtmlFolder() {return gene_html_folder_;}
  bool isKeepUnusedPeaks() {return keep_unused_peaks_;}
  bool isOutputMultipleMass() {return output_multiple_mass_;}
  bool isOutputCsvFeatureFile() {return output_csv_feature_file_;}
  int getThreadNum() {return thread_num_;}
  double getEcscoreCutoff() {return ecscore_cutoff_;}
  bool isSearchPrecWindow() {return search_prec_window_;}
  bool isUseSingleScanNoiseLevel() {return use_single_scan_noise_level_;}

  std::string getMzmlFileName() {return mzml_file_name_;}
  std::string getOutputBaseName() {return output_base_name_;}
  std::string getHtmlDir() {return html_dir_;}
  std::string getMs1JsonDir() {return ms1_json_dir_;}
  std::string getMs2JsonDir() {return ms2_json_dir_;}

  int getFracId() {return frac_id_;}
  bool isFaims() {return is_faims_;}
  double getFaimsVoltage() {return faims_volt_;}
  int getMs1ScanNum() {return ms_1_scan_num_;}
  int getMs2ScanNum() {return ms_2_scan_num_;}
  int getMinScanNum() {return min_scan_num_;}

  void setExeDir(std::string dir) {exe_dir_ = dir;}
  void setResourceDir(std::string dir) {resource_dir_ = dir;}
  void setMissingLevelOne(bool missing) {missing_level_one_ = missing;}
  void setMaxCharge(int charge) {max_charge_ = charge;}
  void setMaxMass(double mass) {max_mass_ = mass;}
  void setMzError(double error) {mz_error_ = error;}
  void setMsOneSnRatio(double ratio) {ms_one_sn_ratio_ = ratio;}
  void setMsTwoSnRatio(double ratio) {ms_two_sn_ratio_ = ratio;}
  void setPrecWindowWidth(double window) { prec_window_ = window;}
  void setUseMsDeconv(bool use) {use_msdeconv_ = use;}
  void setDoFinalFiltering(bool filtering) {do_final_filtering_ = filtering;}
  void setActivation(std::string activation) {activation_ = activation;}
  void setGeneHtmlFolder(bool gene) {gene_html_folder_ = gene;}
  void setKeepUnusedPeaks(bool keep) {keep_unused_peaks_ = keep;}
  void setOutputMultipleMass(bool output) {output_multiple_mass_ = output;}
  void setThreadNum(int num) {thread_num_ = num;}
  void setSearchPrecWindow(bool search) {search_prec_window_ = search;}
  void setUseSingleScanNoiseLevel(bool single_scan_noise)
  {use_single_scan_noise_level_ = single_scan_noise;}
  void setEcscoreCutoff(double cutoff) {ecscore_cutoff_ = cutoff;}
  void setMinScanNum(int min_scan_num) {min_scan_num_ = min_scan_num;}

  void setFracId(int frac_id) {frac_id_ = frac_id;}
  void setMzmlFileNameAndFaims(std::string &mzml_file_name, bool is_faims, double voltage);
  void setMs1ScanNumber(int ms1_scan_num) {ms_1_scan_num_ = ms1_scan_num;}
  void setMs2ScanNumber(int ms2_scan_num) {ms_2_scan_num_ = ms2_scan_num;}

private:
  std::string exe_dir_;
  std::string resource_dir_;
  int max_charge_ = 30;
  double max_mass_ = 70000;
  double prec_window_ = 3.0;
  bool missing_level_one_ = false;
  double mz_error_ = 0.02;
  double ms_one_sn_ratio_ = 3.0;
  double ms_two_sn_ratio_ = 1.0;
  bool keep_unused_peaks_ = false;
  bool use_msdeconv_ = false;
  bool do_final_filtering_ = true;
  int  thread_num_ = 1;
  std::string activation_ = "FILE";
  bool gene_html_folder_ = true;
  double ecscore_cutoff_ = 0.5;
  bool search_prec_window_ = true;
  bool use_single_scan_noise_level_ = false;
  int min_scan_num_ = 3;

  //** Fixed parameter setting **
  // estimate min intensity using the method in Thrash.
  bool estimate_min_inte_ = true;
  bool output_multiple_mass_ = false;
  bool output_match_env_ = false;
  bool output_csv_feature_file_ = true;

  //** information for each run **
  int frac_id_ = -1;
  std::string mzml_file_name_ = "";
  bool is_faims_ = false;
  double faims_volt_ = -1;
  std::string output_base_name_ = "";
  std::string html_dir_ = "";
  std::string ms1_json_dir_ = "";
  std::string ms2_json_dir_ = "";

  int ms_1_scan_num_ = -1;
  int ms_2_scan_num_ = -1;
};

typedef std::shared_ptr<TopdiaPara> TopdiaParaPtr;

}  // namespace toppic

#endif
