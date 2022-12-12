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

#ifndef TOPPIC_TOPFD_DECONV_PROCESS_HPP_
#define TOPPIC_TOPFD_DECONV_PROCESS_HPP_

#include <map>

#include "common/thread/simple_thread_pool.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/env/env_para.hpp"
#include "topfd/common/topfd_para.hpp"
#include "topfd/msreader/raw_ms_group_faime_reader.hpp"
#include "topfd/deconv/deconv_one_sp.hpp"

namespace toppic {
class DeconvProcess {
 public:
  DeconvProcess(){};
  DeconvProcess(TopfdParaPtr topfd_para_ptr, 
                const std::string &spec_file_name, 
                int frac_id, int thread_num);

  void prepareFileFolder(std::string file_num);

  void process();

  void processSp(RawMsGroupFaimeReaderPtr reader_ptr);

  void processSpMissingLevelOne(RawMsGroupFaimeReaderPtr reader_ptr);

  EnvParaPtr getEnvParaPtr() {return env_para_ptr_;}

  std::string updateMsg(MsHeaderPtr header_ptr, int scan, int total_scan_num);

  static int ms1_spec_num_;
  static int ms2_spec_num_;

  int msalign_num_ = 1;//number of ms1/ms2 msalign. it is > 1 if FAIME data.   
  bool isFaims_ = false;//if data is FAIMS data or not

  std::vector<std::pair<double, int>> voltage_vec_;  //{voltage, count of scans with that voltage}
  //store voltages that have appeared in the msgroup so far along with the number of scans with that voltage so far

 private:

  EnvParaPtr env_para_ptr_;
  DpParaPtr dp_para_ptr_;
  TopfdParaPtr topfd_para_ptr_;
  
  std::string spec_file_name_;
  int frac_id_;
  int thread_num_;

  //std::string ms1_env_name_;
  //std::string ms2_env_name_;
  std::string html_dir_;
  std::string ms1_json_dir_;
  std::string ms2_json_dir_; 

};

}

#endif
