//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_DECONV_PROCESS_HPP_
#define PROT_FEATURE_DECONV_PROCESS_HPP_

#include <string>
#include "feature/deconv_para.hpp"
#include "feature/deconv_one_sp.hpp"
#include "feature/feature_mng.hpp"
#include "feature/feature_ms_reader.hpp"

namespace prot {

class DeconvProcess {
 public:
  DeconvProcess(DeconvParaPtr para_ptr) {para_ptr_ = para_ptr;}

  void process();

  void processSp(DeconvOneSpPtr deconv_ptr, FeatureMsReaderPtr reader_ptr, 
                 std::ofstream & ms1_msalign_of, std::ofstream & ms2_msalign_of);

  static void outputParameter(std::ostream &output, DeconvParaPtr para_ptr, const std::string & prefix = "");

 private:
  DeconvParaPtr para_ptr_;

  void copyParameters(FeatureMngPtr mng_ptr);

  std::string updateMsg(MsHeaderPtr header_ptr, int scan, int total_scan_num);
};

}
#endif
