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

#include "prsm/prsm_decoy_filter.hpp"

#include <cstddef>
#include <string>

#include "common/util/file_util.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace toppic {

namespace prsm_decoy_filter {

void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext) {
  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;
  PrsmStrPtrVec prsm_str_ptrs =
      prsm_reader_util::readAllPrsmStrs(input_file_name);
  PrsmStrPtrVec target_ptrs;
  for (size_t i = 0; i < prsm_str_ptrs.size(); i++) {
    if (prsm_str_ptrs[i]->getSeqName().find("DECOY_") != 0) {
      target_ptrs.push_back(prsm_str_ptrs[i]);
    }
  }
  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(target_ptrs);
  writer.close();
}

}  // namespace prsm_decoy_filter

}  // namespace toppic
