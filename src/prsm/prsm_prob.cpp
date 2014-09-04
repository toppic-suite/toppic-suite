#include "base/file_util.hpp"
#include "prsm/prsm_prob.hpp"

namespace prot {

PrsmProb::PrsmProb(const std::string &db_file_name, 
                   const std::string &spec_file_name, 
                   const std::string &in_file_ext,
                   const std::string &out_file_ext, 
                   double K1, double K2, double pref, double inte) {
  input_file_ext_ = in_file_ext;
  output_file_ext_ = out_file_ext;
  spec_file_name_ = spec_file_name;
  db_file_name_ = db_file_name;
  K1_ = K1;
  K2_ = K2;
  pref_ = pref;
  inte_ = inte;
}

void PrsmProb::process() {
  ProteoformPtrVec proteo_ptrs 
      = readFastaToProteoform(db_file_name_, ResidueFactory::getBaseResiduePtrVec());
  LOG_DEBUG("sequences loaded")
  PrsmPtrVec prsm_ptrs 
      = readPrsm(basename(spec_file_name_)+"."+input_file_ext_, proteo_ptrs);
  LOG_DEBUG("prsms loaded")
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int shift_num = prsm_ptrs[i]->getProteoformPtr()->getUnexpectedChangeNum();
    SemiAlignTypePtr type_ptr = prsm_ptrs[i]->getProteoformPtr()->getSemiAlignType();
    ExtremeValuePtr prob_ptr = prsm_ptrs[i]->getProbPtr();
    if (shift_num == 1) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * K1_);
    }
    if (shift_num == 2) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * K2_);
    }
    if (type_ptr == SemiAlignTypeFactory::getPrefixPtr() ||
        type_ptr == SemiAlignTypeFactory::getSuffixPtr() ) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * pref_);
    }

    if (type_ptr == SemiAlignTypeFactory::getInternalPtr()) {
      prob_ptr->setOneProtProb(prob_ptr->getOneProtProb() * inte_);
    }
    
  }
  PrsmWriter all_writer(basename(spec_file_name_)+"."+output_file_ext_);
  all_writer.writeVector(prsm_ptrs);

  //because the prsm_writer ~PrsmWriter changed and the fileclosing is an independant function
  all_writer.close();
}

} /* namespace prot */
