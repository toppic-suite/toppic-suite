/*
 * output_selector.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: xunlikun
 */

#include "prsm/output_selector.hpp"

namespace prot {

OutputSelector::OutputSelector(const std::string &db_file,
                               const std::string &spec_file,
                               const std::string &input_file,
                               const std::string &output_file,
                               const std::string &cutoff_type,
                               double cutoff_value) {
  db_file_= db_file;
  spec_file_ = spec_file;
  input_file_=input_file;
  output_file_ = output_file;
  cutoff_type_ = cutoff_type;
  cutoff_value_ = cutoff_value;
}

OutputSelector::OutputSelector(std::map<std::string,std::string> &arguments,
                               const std::string &input_file,
                               const std::string &output_file) {
  db_file_= arguments["databaseFileName"];
  spec_file_ = arguments["spectrumFileName"];
  input_file_=input_file;
  output_file_ = output_file;
  cutoff_type_ = arguments["cutoffType"];
  cutoff_value_ = atof(arguments["cutoffValue"].c_str());
}

void OutputSelector::process(){
  std::string base_name = basename(spec_file_);
  std::string input_file_name = base_name + "." + input_file_;
  ProteoformPtrVec proteoforms = readFastaToProteoform(db_file_,ResidueFactory::getBaseResiduePtrVec());
  PrSMPtrVec prsms = readPrsm(input_file_name,proteoforms);
  // it's no need to process prsm
  //select
  sort(prsms.begin(),prsms.end(),prsmSpectrumIdUpEvalueUp);
  bool evalue_cutoff = (cutoff_type_ == "EVALUE");

  PrSMPtrVec selected_prsm ;
  ProteoformPtrVec selected_protein;
  int id =0;
  for(unsigned int i=0;i<prsms.size();i++){
    if(evalue_cutoff && prsms[i]->getEValue() <= cutoff_value_){
      prsms[i]->setId(id);
      selected_prsm.push_back(prsms[i]);
      selected_protein.push_back(prsms[i]->getProteoformPtr());
      id++;
    }
    else if(!evalue_cutoff && prsms[i]->getFdr() <= cutoff_value_){
      prsms[i]->setId(id);
      selected_prsm.push_back(prsms[i]);
      selected_protein.push_back(prsms[i]->getProteoformPtr());
      id++;
    }
  }

  //species id

  //output
  std::string output_file_name = base_name +"."+output_file_;
  PrSMWriter writer(output_file_name);
  writer.writeVector(selected_prsm);

  //because the prsm_writer ~PrSMWriter changed and the fileclosing is an independant function
  writer.close();
}
} /* namespace prot */
