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

#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_str_combine.hpp"

namespace toppic {

PrsmStrCombine::PrsmStrCombine(const std::string &spec_file_name, 
                               const std::vector<std::string> &in_file_exts,
                               const std::string &out_file_ext,
                               int top_num):
    spec_file_name_(spec_file_name),
    input_file_exts_(in_file_exts),
    output_file_ext_(out_file_ext),
    top_num_(top_num) {}

PrsmStrCombine::PrsmStrCombine(const std::string &spec_file_name,
                               const std::string &in_file_ext,
                               int in_num,
                               const std::string &out_file_ext,
                               int top_num) {
  output_file_ext_ = out_file_ext;
  spec_file_name_ = spec_file_name;
  top_num_ = top_num;
  for (int i = 0; i < in_num; i ++) {
    std::string ext = in_file_ext + "_" + str_util::toString(i);
    input_file_exts_.push_back(ext);
  }
}

void PrsmStrCombine::process(bool norm) {
  size_t input_num = input_file_exts_.size();
  std::string base_name = file_util::basename(spec_file_name_);
  // open files
  PrsmReaderPtrVec reader_ptrs;
  PrsmStrPtrVec prsm_str_ptrs;
  for (size_t i = 0; i < input_num; i++) {
    std::string input_file_name = base_name + "." + input_file_exts_[i];
    PrsmReaderPtr reader_ptr = std::make_shared<PrsmReader>(input_file_name);
    PrsmStrPtr str_ptr = reader_ptr->readOnePrsmStr();
    reader_ptrs.push_back(reader_ptr);
    prsm_str_ptrs.push_back(str_ptr);
  }
  PrsmXmlWriter writer(base_name + "." + output_file_ext_);

  // combine
  int spec_id = 0;
  bool finish = false;
  while (!finish) {
    finish = true;
    PrsmStrPtrVec cur_str_ptrs;
    for (size_t i = 0; i < input_num; i++) {
      if (prsm_str_ptrs[i] != nullptr) {
        finish = false;
        if (prsm_str_ptrs[i] != nullptr 
            && prsm_str_ptrs[i]->getSpectrumId() < spec_id) {
          LOG_ERROR("Error in the order of reported PrSMs!");
          exit(1);
        }
        while (prsm_str_ptrs[i] != nullptr 
               && prsm_str_ptrs[i]->getSpectrumId() == spec_id) {
          cur_str_ptrs.push_back(prsm_str_ptrs[i]);
          prsm_str_ptrs[i] = reader_ptrs[i]->readOnePrsmStr();
        }
      }
    }

    if (cur_str_ptrs.size() > 0) {
      if (!norm) {
        std::sort(cur_str_ptrs.begin(), cur_str_ptrs.end(), PrsmStr::cmpMatchFragmentDec);
      } else {
        std::sort(cur_str_ptrs.begin(), cur_str_ptrs.end(), PrsmStr::cmpNormMatchFragmentDec);
        auto it = std::unique(cur_str_ptrs.begin(), cur_str_ptrs.end(),
                              [](const PrsmStrPtr & a, const PrsmStrPtr & b) {
                                return a->getSeqName() == b->getSeqName();
                              });
        cur_str_ptrs.erase(it, cur_str_ptrs.end());
      }
      for (unsigned i = 0; i < top_num_; i++) {
        if (i >= cur_str_ptrs.size()) {
          break;
        }
        writer.write(cur_str_ptrs[i]);
      }
    }
    spec_id++;
    LOG_DEBUG("spec id " << spec_id);
  }

  // close files
  for (size_t i = 0; i < input_num; i++) {
    reader_ptrs[i]->close();
  }

  writer.close();
}

} /* namespace toppic */
