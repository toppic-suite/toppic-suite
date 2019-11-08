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


#ifndef TOPPIC_TOPFD_DECONV_PROCESS_HPP_
#define TOPPIC_TOPFD_DECONV_PROCESS_HPP_

#include "ms/spec/msalign_writer.hpp"
#include "topfd/common/topfd_para.hpp"
#include "topfd/msreader/raw_ms_group_reader.hpp"
#include "ms/env/env_para.hpp"
#include "topfd/deconv/deconv_one_sp.hpp"

namespace toppic {

class DeconvProcess {
 public:
  DeconvProcess(TopfdParaPtr topfd_para_ptr, 
                const std::string &spec_file_name, 
                int frac_id);

  void prepareFileFolder();

  void process();

  void processSp(DeconvOneSpPtr deconv_ptr, RawMsGroupReaderPtr reader_ptr, 
                 MsAlignWriterPtr ms1_writer_ptr, MsAlignWriterPtr ms2_writer_ptr);

  void processSpMissingLevelOne(DeconvOneSpPtr deconv_ptr, RawMsGroupReaderPtr reader_ptr, 
                                MsAlignWriterPtr ms2_writer_ptr);

 private:
  EnvParaPtr env_para_ptr_;
  DpParaPtr dp_para_ptr_;
  TopfdParaPtr topfd_para_ptr_;

  //std::string argu_str_;
  //bool missing_level_one_;
  //bool output_match_env_ = false;
  //bool output_json_files_ = true;

  std::string spec_file_name_;
  int frac_id_;

  std::string base_name_;
  std::string ms1_env_name_;
  std::string ms2_env_name_;
  std::string html_dir_;
  std::string ms1_json_dir_;
  std::string ms2_json_dir_; 

  std::string updateMsg(MsHeaderPtr header_ptr, int scan, int total_scan_num);

  void deconvMsOne(RawMsPtr ms_ptr, DeconvOneSpPtr deconv_ptr, 
                   MatchEnvPtrVec &prec_envs, MsAlignWriterPtr ms1_writer_ptr); 

  void deconvMsTwo(RawMsPtr ms_ptr, DeconvOneSpPtr deconv_ptr, 
                   MsAlignWriterPtr ms2_writer_ptr); 
};

}
#endif
