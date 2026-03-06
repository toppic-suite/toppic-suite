//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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
//

#ifndef TOPPIC_TOPFD_DECONV_DECONV_UTIL_HPP_
#define TOPPIC_TOPFD_DECONV_DECONV_UTIL_HPP_

#include <string>

#include "ms/spec/ms_header.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

namespace deconv_util {

std::string updateMsOneMsg(MsHeaderPtr header_ptr, 
                           int scan_cnt, int total_scan_num);

void prepareFileFolder(TopfdParaPtr topfd_para_ptr);

void mergeMs1MsalignFiles(TopfdParaPtr topfd_para_ptr,
                          std::string output_base_name);

}

} // namespace toppic

#endif