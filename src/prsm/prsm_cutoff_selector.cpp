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


#include "base/file_util.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_cutoff_selector.hpp"


namespace prot {

PrsmCutoffSelector::PrsmCutoffSelector(const std::string &db_file_name,
                                       const std::string &spec_file_name,
                                       const std::string &input_file_ext,
                                       const std::string &output_file_ext,
                                       const std::string &cutoff_type,
                                       double cutoff_value): 
    db_file_name_(db_file_name),
    spec_file_name_(spec_file_name),
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext),
    cutoff_type_(cutoff_type),
    cutoff_value_(cutoff_value) {
    }

void PrsmCutoffSelector::process(){
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;

  //PrsmStrPtrVec prsms = PrsmReader::readAllPrsmStrs(input_file_name);
  ModPtrVec fix_mod_list;
  PrsmPtrVec prsms = PrsmReader::readAllPrsms(input_file_name, db_file_name_, fix_mod_list );

  sort(prsms.begin(),prsms.end(),Prsm::cmpSpectrumIdIncPrecursorIdInc);

  bool evalue_cutoff = (cutoff_type_ == "EVALUE");
  bool fdr_cutoff = (cutoff_type_ == "FDR");
  bool form_fdr_cutoff = (cutoff_type_ == "FORMFDR");
  bool frag_cutoff = (cutoff_type_ == "FRAG");

  PrsmPtrVec selected_prsms;
  int id =0; 
  for(size_t i=0; i<prsms.size(); i++){
    if(evalue_cutoff && prsms[i]->getEValue() <= cutoff_value_){
      prsms[i]->setPrsmId(id);
      selected_prsms.push_back(prsms[i]);
      id++;
    }   
    else if(fdr_cutoff && prsms[i]->getFdr() <= cutoff_value_){
      prsms[i]->setPrsmId(id);
      selected_prsms.push_back(prsms[i]);
      id++;
    }   
    else if (form_fdr_cutoff && prsms[i]->getProteoformFdr() <= cutoff_value_) {
      prsms[i]->setPrsmId(id);
      selected_prsms.push_back(prsms[i]);
      id++;
    }
    else if (frag_cutoff && prsms[i]->getMatchFragNum() >= cutoff_value_) {
      prsms[i]->setPrsmId(id);
      ExtremeValuePtr ev_ptr(new ExtremeValue(-prsms[i]->getMatchFragNum(), 1, 1));
      prsms[i]->setExtremeValuePtr(ev_ptr);
      selected_prsms.push_back(prsms[i]);
      id++;
    }   
  }

  //output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(selected_prsms);
  writer.close();
}

} /* namespace prot */
