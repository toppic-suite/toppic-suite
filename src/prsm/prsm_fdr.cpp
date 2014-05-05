/*
 * prsmfdr.cpp
 *
 *  Created on: Mar 20, 2014
 *      Author: xunlikun
 */
#include "base/file_util.hpp"
#include "prsm/prsm_fdr.hpp"


namespace prot {
PrsmFdr::PrsmFdr(std::string db_file_name, std::string spec_file_name,
                 std::string input_ext,std::string output_ext){
  db_file_ = db_file_name;
  spec_file_= spec_file_name;
  input_file_=input_ext;
  output_file_=output_ext;
}
PrsmFdr::PrsmFdr(std::map<std::string,std::string> arguments,std::string input_ext,std::string output_ext){
  db_file_ = arguments["databaseFileName"];
  spec_file_=arguments["spectrumFileName"];
  input_file_=input_ext;
  output_file_=output_ext;
}
void PrsmFdr::process(){
  std::string base_name = basename(spec_file_);
  std::string input_file_name = base_name+"."+input_file_;
  std::string output_file_name = base_name+"."+output_file_;

  ProteoformPtrVec proteoforms = prot::readFastaToProteoform(db_file_,ResidueFactory::getBaseResiduePtrVec());
  PrsmPtrVec prsms = readPrsm(input_file_name,proteoforms);
  PrsmPtrVec target;
  PrsmPtrVec decoy;
  for(unsigned int i=0;i<prsms.size();i++){
    if (prsms[i]->getEValue() == 0.0) {
      LOG_ERROR("prot::PRSMFdr zero E value is reported");
    }
    else {
      std::string prsm_name  = prsms[i]->getProteoformPtr()->getDbResSeqPtr()->getName();
      if(prsm_name.find("DECOY_")==0){
        decoy.push_back(prsms[i]);
      }
      else{
        target.push_back(prsms[i]);
      }
    }
  }
  compute(target,decoy);
  PrsmWriter writer(output_file_name);
  std::sort(target.begin(),target.end(),prsmSpectrumIdUpMatchFragUp);
  writer.writeVector(target);
  writer.close();
}

void PrsmFdr::compute(PrsmPtrVec & target,PrsmPtrVec decoy){
  std::sort(target.begin(),target.end(),prsmEValueUp);
  std::sort(decoy.begin(),decoy.end(),prsmEValueUp);

  /*
  for (unsigned int i=0; i<target.size(); i++){
    LOG_DEBUG("Target " << i << " e value " << target[i]->getEValue());
  }
  for (unsigned int i = 0; i < decoy.size(); i++) {
    LOG_DEBUG("Decoy " << i << " e value " << decoy[i]->getEValue());
  }
  */
  
  for(unsigned int i=0;i<target.size();i++){
    int n_target=i+1;
    double target_evalue = target[i]->getEValue();
    int n_decoy = 0;
    for(unsigned int j=0;j<decoy.size();j++){
      if(decoy[j]->getEValue()<=target_evalue){
        n_decoy++;
      }
      else{
        break;
      }
    }
    //LOG_DEBUG("Target " << i << " e value " << target[i]->getEValue() 
    //          << " target number " << n_target << " decoy number " << n_decoy);
    double fdr = (double)n_decoy/(double)n_target;
    if(fdr>1){
      fdr=1.0;
    }
    target[i]->setFdr(fdr);
  }
}
} /* namespace prot */
