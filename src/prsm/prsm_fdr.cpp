/*
 * prsmfdr.cpp
 *
 *  Created on: Mar 20, 2014
 *      Author: xunlikun
 */
#include "base/file_util.hpp"
#include "prsm/prsm_fdr.hpp"


namespace prot {
PrSMFdr::PrSMFdr(std::string db_file_name, std::string spec_file_name,
                 std::string input_ext,std::string output_ext){
  db_file_ = db_file_name;
  spec_file_= spec_file_name;
  input_file_=input_ext;
  output_file_=output_ext;
}
PrSMFdr::PrSMFdr(std::map<std::string,std::string> arguments,std::string input_ext,std::string output_ext){
  db_file_ = arguments["databaseFileName"];
  spec_file_=arguments["spectrumFileName"];
  input_file_=input_ext;
  output_file_=output_ext;
}
void PrSMFdr::process(){
  std::string base_name = basename(spec_file_);
  std::string input_file_name = base_name+"."+input_file_;
  std::string output_file_name = base_name+"."+output_file_;

  ProteoformPtrVec proteoforms = prot::readFastaToProteoform(db_file_,ResidueFactory::getBaseResiduePtrVec());
  PrSMPtrVec prsms = readPrsm(input_file_name,proteoforms);
  PrSMPtrVec target;
  PrSMPtrVec decoy;
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
  PrSMWriter writer(output_file_name);
  std::sort(target.begin(),target.end(),prsmSpectrumIdUpMatchFragUp);
  writer.writeVector(target);
  writer.close();
}

void PrSMFdr::compute(PrSMPtrVec & target,PrSMPtrVec decoy){
  std::sort(target.begin(),target.end(),prsmEValueUp);
  std::sort(decoy.begin(),decoy.end(),prsmEValueUp);

  for (unsigned int i=0; i<target.size(); i++){
    std::cout << "Target " << i << " e value " << target[i]->getEValue() << std::endl;
  }

  for (unsigned int i = 0; i < decoy.size(); i++) {
    std::cout << "Decoy " << i << " e value " << decoy[i]->getEValue() << std::endl;
  }
  
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
    std::cout << "Target " << i << " e value " << target[i]->getEValue() << " target number " << n_target << " decoy number " << n_decoy << std::endl;
    double fdr = (double)n_decoy/(double)n_target;
    if(fdr>1){
      fdr=1.0;
    }
    target[i]->setFdr(fdr);
  }
}
} /* namespace prot */
