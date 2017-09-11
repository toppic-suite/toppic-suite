// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <string>
#include <algorithm>

#include "base/file_util.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_form_filter.hpp"


namespace prot {

void PrsmFormFilter::process() {
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;

  // PrsmStrPtrVec prsms = PrsmReader::readAllPrsmStrs(input_file_name);
  ModPtrVec fix_mod_list;
  PrsmPtrVec prsms = PrsmReader::readAllPrsms(input_file_name, db_file_name_, fix_mod_list);

  sort(prsms.begin(), prsms.end(), Prsm::cmpEValueInc);

  PrsmPtrVec selected_prsms;

  PrsmPtrVec selected_forms;

  for (size_t i = 0; i < prsms.size(); i++) {
    // std::cout << "prsm " << i << std::endl;
    bool found = false;
    for (size_t j = 0; j < selected_forms.size(); j++) {
      if (selected_forms[j]->getProteoformPtr()->getSpeciesId()
          == prsms[i]->getProteoformPtr()->getSpeciesId()) {
        found = true;
        break;
      }
    }
    if (found) {
      selected_prsms.push_back(prsms[i]);
    } else {
      bool keep = true;
      std::string form = prsms[i]->getProteoformPtr()->getProteinMatchSeq();
      for (size_t j = 0; j < selected_forms.size(); j++) {
        if (selected_forms[j]->getProteoformPtr()->getSpeciesId() == prsms[i]->getProteoformPtr()->getSpeciesId()) {
          // std::cout << "scan " << prsms[i]->getSpectrumScan() << " removed by scan " << selected_forms[j]->getSpectrumScan() << std::endl;
          keep = false;
          break;
        }
        if (selected_forms[j]->getProteoformPtr()->getProteinMatchSeq() == form) {
          // std::cout << "scan " << prsms[i]->getSpectrumScan() << " removed by scan " << selected_forms[j]->getSpectrumScan() << std::endl;
          keep = false;
          break;
        }
      }
      if (keep) {
        selected_forms.push_back(prsms[i]);
        selected_prsms.push_back(prsms[i]);
      }
    }
  }

  // output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  sort(selected_prsms.begin(), selected_prsms.end(), Prsm::cmpSpectrumIdIncPrecursorIdInc);
  writer.writeVector(selected_prsms);
  writer.close();

  output_file_name = base_name + "." + output_file_ext_2_;
  PrsmXmlWriter writer2(output_file_name);
  sort(selected_forms.begin(), selected_forms.end(), Prsm::cmpSpectrumIdIncPrecursorIdInc);
  writer2.writeVector(selected_forms);
  writer2.close();
}

} /* namespace prot */
