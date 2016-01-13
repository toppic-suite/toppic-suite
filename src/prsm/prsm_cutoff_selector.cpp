#include "base/file_util.hpp"
#include "prsm/prsm_str.hpp"
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

  PrsmStrPtrVec prsms = PrsmReader::readAllPrsmStrs(input_file_name);

  //select
  sort(prsms.begin(),prsms.end(),PrsmStr::cmpSpectrumIdInc);
  bool evalue_cutoff = (cutoff_type_ == "EVALUE");

  PrsmStrPtrVec selected_prsms;
  int id =0;
  for(size_t i=0; i<prsms.size(); i++){
    if(evalue_cutoff && prsms[i]->getEValue() <= cutoff_value_){
      prsms[i]->setId(id);
      selected_prsms.push_back(prsms[i]);
      id++;
    }
    else if(!evalue_cutoff && prsms[i]->getFdr() <= cutoff_value_){
      prsms[i]->setId(id);
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
