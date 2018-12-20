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


#ifndef PROT_PRSM_SIMPLE_PRSM_STR_COMBINE_HPP_
#define PROT_PRSM_SIMPLE_PRSM_STR_COMBINE_HPP_

#include <vector>
#include <string>
#include <map>

#include "seq/proteoform.hpp"
#include "seq/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"

namespace toppic {

class SimplePrsmStrCombine {
 public:
  SimplePrsmStrCombine(const std::string &spec_file_name,
                       const std::vector<std::string> &in_file_exts,
                       const std::string &out_file_ext,
                       int top_num):
      spec_file_name_(spec_file_name),
      input_file_exts_(in_file_exts),
      output_file_ext_(out_file_ext),
      top_num_(top_num) {}

  SimplePrsmStrCombine(const std::string &spec_file_name,
                       const std::string &in_file_ext,
                       int in_num,
                       const std::string &out_file_ext,
                       int top_num):
      spec_file_name_(spec_file_name),
      output_file_ext_(out_file_ext),
      top_num_(top_num) {
        for (int i = 0; i < in_num; i ++) {
          std::string ext = in_file_ext + "_" + str_util::toString(i);
          input_file_exts_.push_back(ext);
        }
      }

  void process();

 private:
  std::string spec_file_name_;
  std::vector<std::string> input_file_exts_;
  std::string output_file_ext_;
  int top_num_;
};

typedef std::shared_ptr<SimplePrsmStrCombine> SimplePrsmStrCombinePtr;
} /* namespace toppic */

#endif
