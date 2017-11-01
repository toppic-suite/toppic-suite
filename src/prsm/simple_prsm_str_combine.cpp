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

#include <set>
#include <string>
#include <algorithm>

#include "base/file_util.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace prot {

void SimplePrsmStrCombine::process() {
  size_t input_num = input_file_exts_.size();
  std::string base_name = FileUtil::basename(spec_file_name_);
  // open files
  SimplePrsmReaderPtrVec reader_ptrs;
  SimplePrsmStrPtrVec prsm_str_ptrs;
  for (size_t i = 0; i < input_num; i++) {
    std::string input_file_name = base_name + "." + input_file_exts_[i];
    SimplePrsmReaderPtr reader_ptr
        = std::make_shared<SimplePrsmReader>(input_file_name);
    LOG_DEBUG("input file name " << input_file_name);
    SimplePrsmStrPtr str_ptr = reader_ptr->readOnePrsmStr();
    reader_ptrs.push_back(reader_ptr);
    prsm_str_ptrs.push_back(str_ptr);
  }
  SimplePrsmXmlWriter writer(base_name + "." + output_file_ext_);

  // combine
  int spec_id = 0;
  bool finish = false;
  while (!finish) {
    // LOG_DEBUG("spec id " << spec_id << " input num " << input_num);
    finish = true;
    SimplePrsmStrPtrVec cur_str_ptrs;
    for (size_t i = 0; i < input_num; i++) {
      if (prsm_str_ptrs[i] != nullptr) {
        finish = false;
        while (prsm_str_ptrs[i]!= nullptr &&
               prsm_str_ptrs[i]->getSpectrumId() == spec_id) {
          cur_str_ptrs.push_back(prsm_str_ptrs[i]);
          prsm_str_ptrs[i] = reader_ptrs[i]->readOnePrsmStr();
        }
      }
    }
    // LOG_DEBUG("finish " << finish);

    if (cur_str_ptrs.size() > 0) {
      std::sort(cur_str_ptrs.begin(), cur_str_ptrs.end(), SimplePrsmStr::cmpScoreDec);
      int count = 0;
      std::set<std::string> name_set;
      for (size_t i = 0; i < cur_str_ptrs.size(); i++) {
        if (count >= top_num_) break;
        if (name_set.find(cur_str_ptrs[i]->getSeqName()) == name_set.end()) {
          count++;
          name_set.insert(cur_str_ptrs[i]->getSeqName());
          writer.write(cur_str_ptrs[i]);
        }
      }
    }
    spec_id++;
  }

  // close files
  for (size_t i = 0; i < input_num; i++) {
    reader_ptrs[i]->close();
  }

  writer.close();
}

} /* namespace prot */
