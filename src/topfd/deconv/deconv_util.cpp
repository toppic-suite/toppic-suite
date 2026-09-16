// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include "topfd/deconv/deconv_util.hpp"

#include <iostream>

#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/spec/msalign_thread_merge.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

namespace deconv_util {

std::string updateMsOneMsg(const MsHeaderPtr& header_ptr, int scan_cnt,
                           int total_scan_num) {
  std::string percentage = std::to_string(scan_cnt * 100 / total_scan_num);
  std::string msg = "Processing MS1 spectrum scan " +
                    std::to_string(header_ptr->getFirstScanNum()) + " ...";
  while (msg.length() < 40) {
    msg += " ";
  }
  msg = msg + percentage + "% finished.";
  return msg;
}

void mergeMs1MsalignFiles(const TopfdParaPtr& topfd_para_ptr,
                          std::string output_base_name) {
  // Merge files
  std::string file_name_ext = "ms1.msalign";
  std::string para_str = topfd_para_ptr->getParaStr("#", "\t");
  MsalignThreadMergePtr ms1_merge_ptr = std::make_shared<MsalignThreadMerge>(
      file_name_ext, topfd_para_ptr->getThreadNum(), file_name_ext,
      output_base_name, para_str);
  ms1_merge_ptr->process();

  // remove temporary files
  std::string ms1_prefix =
      file_util::absoluteName(output_base_name) + "_ms1.msalign_";
  file_util::cleanPrefix(output_base_name, ms1_prefix);
  std::cout << std::endl;
}

}  // namespace deconv_util

};  // namespace toppic
