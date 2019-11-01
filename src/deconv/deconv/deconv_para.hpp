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


#ifndef TOPPIC_DECONV_DECONV_DECONV_PARA_HPP_
#define TOPPIC_DECONV_DECONV_DECONV_PARA_HPP_

#include <memory>
#include <string>
#include <map>

namespace toppic {

class DeconvPara {
 public:
  DeconvPara(std::map<std::string, std::string> &arguments, 
             const std::string &argument_str); 

  std::string resource_dir_;

  bool refine_prec_mass_;

  bool missing_level_one_;

  int max_charge_;

  double max_mass_;

  double tolerance_;

  double ms_two_sn_ratio_;

  double ms_one_sn_ratio_;

  bool keep_unused_peaks_;

  bool output_multiple_mass_ = false; 

  double prec_window_;
  
  bool do_final_filtering_ = true;

  bool output_match_env_ = false;

  bool output_json_files_ = true;

  std::string argument_str_;
};

typedef std::shared_ptr<DeconvPara> DeconvParaPtr;

}
#endif
