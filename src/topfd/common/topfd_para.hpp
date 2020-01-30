//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

  std::string getParaStr(const std::string &prefix);

  std::string resource_dir_;
  bool refine_prec_mass_ = true;
  bool missing_level_one_ = false;
  int max_charge_ = 30;
  double max_mass_ = 100000;
  double mz_error_ = 0.02;
  double ms_one_sn_ratio_ = 3.0;
  double ms_two_sn_ratio_ = 1.0;
  double prec_window_ = 3.0;
  bool keep_unused_peaks_ = false;
  bool output_multiple_mass_ = false;
  bool do_final_filtering_ = true;
  bool output_match_env_ = false;
  bool output_json_files_ = true;
  bool merge_files_ = false;
  std::string thread_number = "1";
  std::string merged_file_name_ = "";
};

typedef std::shared_ptr<TopfdPara> TopfdParaPtr;

}  // namespace toppic

#endif 
