//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_TOPFD_COMMON_TOPFD_PARA_HPP_
#define TOPPIC_TOPFD_COMMON_TOPFD_PARA_HPP_

#include <memory>
#include <string>

namespace toppic {

class TopfdPara {
 public:
  TopfdPara() {};
  
  std::string getParaStr(const std::string &prefix,
		         const std::string &sep);

  std::string getExeDir() {return exe_dir_;}
  std::string getResourceDir() {return resource_dir_;}
  bool isMissingLevelOne() {return missing_level_one_;}
  int getMaxCharge() {return max_charge_;}
  double getMaxMass() {return max_mass_;}
  double getMzError() {return mz_error_;}
  double getMsOneSnRatio() {return ms_one_sn_ratio_;}
  double getMsTwoSnRatio() {return ms_two_sn_ratio_;}
  double getPrecWindow() {return prec_window_;}
  bool isUseEnvCnn() {return use_env_cnn_;}
  bool isDoFinalFiltering() {return do_final_filtering_;}
  std::string getActivation() {return activation_;}
  bool isGeneHtmlFolder() {return gene_html_folder_;}
  bool isKeepUnusedPeaks() {return keep_unused_peaks_;}
  bool isOutputMultipleMass() {return output_multiple_mass_;}
  int getThreadNum() {return thread_num_;}

  void setExeDir(std::string dir) {exe_dir_ = dir;}
  void setResourceDir(std::string dir) {resource_dir_ = dir;}
  void setMissingLevelOne(bool missing) {missing_level_one_ = missing;}
  void setMaxCharge(int charge) {max_charge_ = charge;}
  void setMaxMass(double mass) {max_mass_ = mass;}
  void setMzError(double error) {mz_error_ = error;}
  void setMsOneSnRatio(double ratio) {ms_one_sn_ratio_ = ratio;}
  void setMsTwoSnRatio(double ratio) {ms_two_sn_ratio_ = ratio;}
  void setPrecWindow(double window) {prec_window_ = window;}
  void setUseEnvCnn(bool use) {use_env_cnn_ = use;}
  void setDoFinalFiltering(bool filtering) {do_final_filtering_ = filtering;}
  void setActivation(std::string activation) {activation_ = activation;}
  void setGeneHtmlFolder(bool gene) {gene_html_folder_ = gene;}
  void setKeepUnusedPeaks(bool keep) {keep_unused_peaks_ = keep;}
  void setOutputMultipleMass(bool output) {output_multiple_mass_ = output;}
  void setThreadNum(int num) {thread_num_ = num;}

  void setMs1ScanNumber(int ms1_scan_num) {ms_1_scan_num_ = ms1_scan_num;}
  void setMs2ScanNumber(int ms2_scan_num) {ms_2_scan_num_ = ms2_scan_num;}

 private:
  std::string exe_dir_;
  std::string resource_dir_;
  bool refine_prec_mass_ = true;
  bool missing_level_one_ = false;
  int max_charge_ = 30;
  double max_mass_ = 70000;
  double mz_error_ = 0.02;
  double ms_one_sn_ratio_ = 3.0;
  double ms_two_sn_ratio_ = 1.0;
  double prec_window_ = 3.0;
  bool keep_unused_peaks_ = false;
  bool use_env_cnn_ = false;
  bool output_multiple_mass_ = false;
  bool do_final_filtering_ = true;
  bool output_match_env_ = false;
  int  thread_num_ = 1;
  std::string activation_ = "FILE";
  bool gene_html_folder_ = true;
  int ms_1_scan_num_ = -1;
  int ms_2_scan_num_ = -1;
};

typedef std::shared_ptr<TopfdPara> TopfdParaPtr;

}  // namespace toppic

#endif 
