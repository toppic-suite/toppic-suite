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

#include <algorithm>

#include "common/util/file_util.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_cutoff_selector.hpp"

namespace toppic {

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
    cutoff_value_(cutoff_value) {}

void PrsmCutoffSelector::process(){
  std::string base_name = file_util::basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;

  ModPtrVec fix_mod_list;
  PrsmPtrVec prsms = PrsmReader::readAllPrsms(input_file_name, db_file_name_, fix_mod_list );

  std::sort(prsms.begin(), prsms.end(), Prsm::cmpSpectrumIdIncPrecursorIdInc);

  bool evalue_cutoff = (cutoff_type_ == "EVALUE");
  bool fdr_cutoff = (cutoff_type_ == "FDR");
  bool form_fdr_cutoff = (cutoff_type_ == "FORMFDR");
  bool frag_cutoff = (cutoff_type_ == "FRAG");

  PrsmPtrVec selected_prsms;
  int id = 0; 
  for(size_t i = 0; i<prsms.size(); i++){
    if(evalue_cutoff && prsms[i]->getEValue() <= cutoff_value_){
      prsms[i]->setPrsmId(id);
      selected_prsms.push_back(prsms[i]);
      id++;
    } else if(fdr_cutoff && prsms[i]->getFdr() <= cutoff_value_){
      prsms[i]->setPrsmId(id);
      selected_prsms.push_back(prsms[i]);
      id++;
    } else if (form_fdr_cutoff && prsms[i]->getFdr() <= cutoff_value_
               && prsms[i]->getProteoformFdr() <= cutoff_value_) {
      prsms[i]->setPrsmId(id);
      selected_prsms.push_back(prsms[i]);
      id++;
    } else if (frag_cutoff && prsms[i]->getMatchFragNum() >= cutoff_value_) {
      prsms[i]->setPrsmId(id);
      ExtremeValuePtr ev_ptr = std::make_shared<ExtremeValue>(-prsms[i]->getMatchFragNum(), 1, 1);
      prsms[i]->setExtremeValuePtr(ev_ptr);
      selected_prsms.push_back(prsms[i]);
      id++;
    }   
  }

  // output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(selected_prsms);
  writer.close();
}

} /* namespace toppic */
