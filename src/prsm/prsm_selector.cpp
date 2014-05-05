/*
 * prsm_selector.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: xunlikun
 */

#include "base/file_util.hpp"
#include "prsm/prsm_selector.hpp"

namespace prot {

PrsmSelector::PrsmSelector(std::string db_file,std::string spec_file,std::string in_file,std::string out_file,int n_top){
  spec_file_ = spec_file;
  db_file_ = db_file;
  input_file_ = in_file;
  output_file_ = out_file;
  n_top_ = n_top;
}

PrsmSelector::PrsmSelector(std::map<std::string,std::string> arguments,std::string in_file,std::string out_file,int n_top){
  spec_file_ = arguments["spectrumFileName"];
  db_file_ = arguments["databaseFileName"];
  input_file_ = in_file;
  output_file_ = out_file;
  n_top_ = n_top;
}
bool PrsmSelector::findPrsm(PrsmPtrVec result,PrsmPtr prsm){
  for(unsigned int i=0;i< result.size();i++){
    if(result[i]->getProteoformPtr()->getDbResSeqPtr()->getId()==
        prsm->getProteoformPtr()->getDbResSeqPtr()->getId()){
      return true;
    }
  }
  return false;
}
PrsmPtrVec PrsmSelector::getTopPrsms(PrsmPtrVec selected_prsm,int n_top){
  std::sort(selected_prsm.begin(),selected_prsm.end(),prsmEValueUp);
  int size = selected_prsm.size();
  int max = size > n_top? n_top:size;
  PrsmPtrVec result;
  for(int i=0;i<max;i++){
    if(!findPrsm(result,selected_prsm[i])){
      result.push_back(selected_prsm[i]);
    }
  }
  return result;
}
void PrsmSelector::process(){
  std::string base_name = basename(spec_file_);
  std::string input_file_name = base_name+"."+input_file_;
  ProteoformPtrVec proteoforms = prot::readFastaToProteoform(db_file_,ResidueFactory::getBaseResiduePtrVec());
  PrsmPtrVec prsms = readPrsm(input_file_name,proteoforms);
  std::string output_file_name = base_name+"."+output_file_;
//  sort(prsms.begin(),prsms.end(),prsmSpectrumIdUpMatchFragUp);
  int max_id = prsms[prsms.size()-1]->getSpectrumId();
  PrsmWriter writer(output_file_name);
  for(int i=0;i<= max_id;i++){
    PrsmPtrVec selected_prsms;
    //select prsm by spectrum id
    for(unsigned int j=0;j<prsms.size();j++){
      if(prsms[j]->getSpectrumId()==i){
        selected_prsms.push_back(prsms[j]);
      }
    }
    PrsmPtrVec result = getTopPrsms(selected_prsms,n_top_);
    writer.writeVector(result);
  }
  //because the prsm_writer ~PrsmWriter changed and the fileclosing is an independant function
  writer.close();
}

} /* namespace prot */
