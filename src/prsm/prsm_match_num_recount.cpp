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

#include "prsm/prsm_match_num_recount.hpp"

#include <cstddef>
#include <string>

#include "common/util/file_util.hpp"
#include "para/prsm_para.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace toppic {

namespace prsm_match_num_recount {

void process(const PrsmParaPtr& prsm_para_ptr,
             const std::string& input_file_ext,
             const std::string& output_file_ext) {
  // attaches the full deconvoluted and refined spectra to each PrSM
  PrsmPtrVec prsm_ptrs =
      prsm_reader_util::readPrsmsWithSpectra(prsm_para_ptr, input_file_ext);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    prsm_ptrs[i]->updateMatchNum(sp_para_ptr);
  }
  std::string output_file_name =
      file_util::basename(prsm_para_ptr->getSpectrumFileName()) + "." +
      output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}

}  // namespace prsm_match_num_recount

}  // namespace toppic
