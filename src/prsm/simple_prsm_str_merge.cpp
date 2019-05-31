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

#include <set>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_str_merge.hpp"

namespace toppic {

SimplePrsmStrMerge::SimplePrsmStrMerge(const std::string &spec_file_name,
                                       const std::vector<std::string> &in_file_exts,
                                       const std::string &out_file_ext,
                                       int top_num):
    spec_file_name_(spec_file_name),
    input_file_exts_(in_file_exts),
    output_file_ext_(out_file_ext),
    top_num_(top_num) {}

SimplePrsmStrMerge::SimplePrsmStrMerge(const std::string &spec_file_name,
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

SimplePrsmStrMerge::SimplePrsmStrMerge(const std::string &spec_file_name,
                                       const std::string &in_file_pref,
                                       const std::string &in_file_suff,
                                       int in_num,
                                       const std::string &out_file_ext,
                                       int top_num):
    spec_file_name_(spec_file_name),
    output_file_ext_(out_file_ext),
    top_num_(top_num) {
      for (int i = 0; i < in_num; i ++) {
        std::string ext = in_file_pref + "_" + str_util::toString(i) + "_" + in_file_suff;
        input_file_exts_.push_back(ext);
      }
    }

void SimplePrsmStrMerge::process() {
  size_t input_num = input_file_exts_.size();
  std::string base_name = file_util::basename(spec_file_name_);
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

void SimplePrsmStrMerge::mergeBlockResults(std::string &sp_file_name, 
                                           std::string &input_pref,
                                           int block_num, 
                                           int comp_num, 
                                           int pref_suff_num,
                                           int inte_num) {

  std::string complete = ProteoformType::COMPLETE->getName();
  std::string prefix = ProteoformType::PREFIX->getName();
  std::string suffix = ProteoformType::SUFFIX->getName();
  std::string internal = ProteoformType::INTERNAL->getName();

  SimplePrsmStrMerge comp_merge(sp_file_name, 
                                input_pref,
                                complete,
                                block_num, 
                                input_pref + "_" + complete,
                                comp_num);
  comp_merge.process();

  SimplePrsmStrMerge pref_merge(sp_file_name, 
                                input_pref,
                                prefix,
                                block_num, 
                                input_pref + "_" + prefix,
                                pref_suff_num);
  pref_merge.process();

  SimplePrsmStrMerge suff_merge(sp_file_name, 
                                input_pref,
                                suffix,
                                block_num, 
                                input_pref + "_" + suffix,
                                pref_suff_num);
  suff_merge.process();

  SimplePrsmStrMerge internal_merge(sp_file_name, 
                                    input_pref, 
                                    internal,
                                    block_num, 
                                    input_pref + "_" + internal,
                                    inte_num);
  internal_merge.process();

  // remove tempory files
  for (int i = 0; i < block_num; i++) {
    std::string file_pref = input_pref + "_" + std::to_string(i) + "_";
    file_util::cleanTempFiles(sp_file_name, file_pref);
  }
}

} /* namespace toppic */
