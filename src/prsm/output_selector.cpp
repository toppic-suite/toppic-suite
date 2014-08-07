#include "base/file_util.hpp"
#include "prsm/output_selector.hpp"

namespace prot {

OutputSelector::OutputSelector(const std::string &db_file_name,
                               const std::string &spec_file_name,
                               const std::string &input_file_ext,
                               const std::string &output_file_ext,
                               const std::string &cutoff_type,
                               double cutoff_value) {
  db_file_name_= db_file_name;
  spec_file_name_ = spec_file_name;
  cutoff_type_ = cutoff_type;
  cutoff_value_ = cutoff_value;
  input_file_ext_= input_file_ext;
  output_file_ext_ = output_file_ext;
}

void OutputSelector::process(){
  std::string base_name = basename(spec_file_name_);
  std::string input_file_name = base_name + "." + input_file_ext_;
  ProteoformPtrVec proteoforms 
      = readFastaToProteoform(db_file_name_,ResidueFactory::getBaseResiduePtrVec());
  PrsmPtrVec prsms = readPrsm(input_file_name,proteoforms);
  // it's no need to process prsm
  //select
  sort(prsms.begin(),prsms.end(),prsmSpectrumIdUpEvalueUp);
  bool evalue_cutoff = (cutoff_type_ == "EVALUE");

  PrsmPtrVec selected_prsms;
  ProteoformPtrVec selected_prots;
  int id =0;
  for(size_t i=0; i<prsms.size(); i++){
    if(evalue_cutoff && prsms[i]->getEValue() <= cutoff_value_){
      prsms[i]->setId(id);
      selected_prsms.push_back(prsms[i]);
      selected_prots.push_back(prsms[i]->getProteoformPtr());
      id++;
    }
    else if(!evalue_cutoff && prsms[i]->getFdr() <= cutoff_value_){
      prsms[i]->setId(id);
      selected_prsms.push_back(prsms[i]);
      selected_prots.push_back(prsms[i]->getProteoformPtr());
      id++;
    }
  }

  //output
  std::string output_file_name = base_name + "." + output_file_ext_;
  PrsmWriter writer(output_file_name);
  writer.writeVector(selected_prsms);

  //because the prsm_writer ~PrsmWriter changed and the fileclosing is an independant function
  writer.close();
}
} /* namespace prot */
