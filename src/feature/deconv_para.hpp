//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_DECONV_PARA_HPP_
#define PROT_FEATURE_DECONV_PARA_HPP_

#include <string>
#include <memory>
#include <map>
#include <vector>

namespace prot {

class DeconvPara {
 public:
  DeconvPara(std::map<std::string, std::string> &arguments);

  void setDataFileName(const std::string & file_name) {data_file_name_ = file_name;}

  std::string getDataFileName() {return data_file_name_;}

  std::string data_file_name_;
  std::string exec_dir_;

  bool refine_prec_mass_;
  bool missing_level_one_;

  int max_charge_;
  double max_mass_;
  double tolerance_;
  double sn_ratio_;
  bool keep_unused_peaks_;

  bool output_multiple_mass_ = false; 

  double prec_window_;
};

typedef std::shared_ptr<DeconvPara> DeconvParaPtr;

}
#endif
