/*
 * prsm_combine.cpp
 *
 *  Created on: Feb 18, 2014
 *      Author: xunlikun
 */

#include "prsm/prsm_combine.hpp"

namespace prot {

PrSMCombine::PrSMCombine(std::string db_file,
                         std::string spec_file,
                         std::vector<std::string> &in_file_exts,
                         std::string out_file) {
  input_file_exts_ = in_file_exts;
  output_files_ = out_file;
  spec_file_ = spec_file;
  db_file_ = db_file;
}

PrSMCombine::~PrSMCombine() {
}

void PrSMCombine::process() {
  ProteoformPtrVec proteoforms = prot::readFastaToProteoform(db_file_,ResidueFactory::getBaseResiduePtrVec());
  PrSMPtrVec prsms;
  for(unsigned int i=0;i<input_file_exts_.size();i++){
    PrSMPtrVec temps = prot::readPrsm(basename(spec_file_)+"."+input_file_exts_[i],proteoforms);
    for(unsigned int j=0;j<temps.size();j++){
//      std::cout<<temps[j]->getSpectrumId()<<std::endl;
      prsms.push_back(temps[j]);
    }
  }
  std::sort(prsms.begin(),prsms.end(),prsm_spectrum);
  PrSMWriter all_writer(basename(spec_file_)+"."+output_files_);
  all_writer.writeVector(prsms);
}

} /* namespace prot */
