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

#include <set>
#include <algorithm>
#include <iostream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/spec/msalign_thread_merge.hpp"

namespace toppic {

MsalignThreadMerge::MsalignThreadMerge(const std::string &spec_file_name,
                                       const std::string &in_file_ext,
                                       int in_num,
                                       const std::string &out_file_ext, const std::string &para_str):
    spec_file_name_(spec_file_name),
    para_str_(para_str),
    output_file_ext_(out_file_ext) {
      for (int i = 0; i < in_num; i ++) {
        std::string ext = in_file_ext + "_" + str_util::toString(i);
        input_file_exts_.push_back(ext);
      }
    }

inline int getCurMsIndex(DeconvMsPtrVec &ms_ptrs) {
  int index = -1;
  int scan_num = std::numeric_limits<int>::max();
  for (size_t i = 0; i < ms_ptrs.size(); i++) {
    if (ms_ptrs[i] != nullptr) {
      if (ms_ptrs[i]->getMsHeaderPtr()->getFirstScanNum() < scan_num) {
        scan_num = ms_ptrs[i]->getMsHeaderPtr()->getFirstScanNum();
        index = i;
      }
    }
  }
  return index;
}

void MsalignThreadMerge::process() {
  size_t input_num = input_file_exts_.size();
  std::string base_name = file_util::basename(spec_file_name_);
  // open files
  MsAlignReaderPtrVec reader_ptrs;
  DeconvMsPtrVec ms_ptrs;
  for (size_t i = 0; i < input_num; i++) {
    std::string input_file_name = base_name + "_" + input_file_exts_[i];
    MsAlignReaderPtr reader_ptr
        = std::make_shared<MsAlignReader>(input_file_name);
    LOG_DEBUG("input file name " << input_file_name);
    DeconvMsPtr ms_ptr = reader_ptr->getNextMs();
    reader_ptrs.push_back(reader_ptr);
    ms_ptrs.push_back(ms_ptr);
  }
  std::string output_filename = base_name + "_" + output_file_ext_;
  MsAlignWriterPtr writer = std::make_shared<MsAlignWriter>(output_filename); 
  
  writer->writePara(para_str_);

  // combine
  int spec_id = 0;
  bool finish = true;
  
  while (finish) {
    // LOG_DEBUG("spec id " << spec_id << " input num " << input_num);
    int cur_ms_idx = getCurMsIndex(ms_ptrs);

    if (cur_ms_idx < 0) {
      break;
    }
    else { 
      DeconvMsPtr cur_ms_ptr = ms_ptrs[cur_ms_idx];
      cur_ms_ptr->getMsHeaderPtr()->setId(spec_id);
      spec_id++;
      writer->write(cur_ms_ptr);
      ms_ptrs[cur_ms_idx] = reader_ptrs[cur_ms_idx]->getNextMs();
    }
  }

  // close files
  for (size_t i = 0; i < input_num; i++) {
    reader_ptrs[i]->close();
  }
  writer->close();
}

} /* namespace toppic */
