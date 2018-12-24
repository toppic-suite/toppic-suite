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


#ifndef TOPPIC_DECONV_DECONV_DECONV_PROCESS_HPP_
#define TOPPIC_DECONV_DECONV_DECONV_PROCESS_HPP_

#include "deconv/msreader/raw_ms_reader.hpp"
#include "deconv/env/env_para.hpp"
#include "deconv/deconv/deconv_para.hpp"
#include "deconv/deconv/deconv_one_sp.hpp"

namespace toppic {

class DeconvProcess {
 public:
  DeconvProcess(DeconvParaPtr para_ptr) {para_ptr_ = para_ptr;}

  void process();

  void processSp(DeconvOneSpPtr deconv_ptr, RawMsReaderPtr reader_ptr, 
                 std::ofstream & ms1_msalign_of, std::ofstream & ms2_msalign_of);

  static std::string getParameterStr(DeconvParaPtr para_ptr, const std::string & prefix = "");

 private:
  DeconvParaPtr para_ptr_;

  void copyParameters(EnvParaPtr env_para_ptr);

  std::string updateMsg(MsHeaderPtr header_ptr, int scan, int total_scan_num);
};

}
#endif
