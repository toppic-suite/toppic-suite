#include <sstream>

#include "base/proteoform.hpp"

#include <log4cxx/logger.h>

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Proteofom"));

Proteoform::Proteoform(ProteoformPtr ori_form_ptr, std::string name, ResSeqPtr res_seq_ptr, 
             int start_pos, int end_pos, ChangePtrVec change_list) {
  ori_form_ptr_ = ori_form_ptr;
  name_ = name;
  residue_seq_ptr_ = res_seq_ptr;
  start_pos_ = start_pos;
  end_pos_ = end_pos;
  LOG4CXX_TRACE(logger, "start bp spec generation");
  bp_spec_ptr_ = BpSpecPtr(new BpSpec(res_seq_ptr));
  change_list_ = change_list;
}


ProteoformPtr Proteoform::getProtModProteoform(ProteoformPtr ori_form_ptr,
                                               ResiduePtrVec &residue_list, 
                                               ProtModPtr prot_mod_ptr) {
  TruncPtr trunc_ptr = prot_mod_ptr->getTruncPtr();
  int len = trunc_ptr->getTruncLen();
  AcidPtrVec trunc_acids = trunc_ptr->getAcidPtrVec();
  //check if trunc acids match N-terminal acids of the protein 
  if (len >= residue_seq_ptr_->getLen()) {
    return ProteoformPtr(nullptr);
  }
  bool match = true;
  for (int i = 0; i < len; i++) {
    AcidPtr acid = residue_seq_ptr_->getResiduePtr(i)->getAcidPtr();
    if (acid.get() != trunc_acids[i].get()) {
      match = false;
      break;
    }
  }
  if (!match) {
    return ProteoformPtr(nullptr);
  }
  ResiduePtrVec residues;
  ResiduePtr residue = residue_seq_ptr_->getResiduePtr(len);
  PtmPtr ori_ptm = residue->getPtmPtr();
  PtmPtr ptm = prot_mod_ptr->getPtmPtr();
  ChangePtrVec change_list;
  if (ptm->isEmpty() || !ori_ptm->isEmpty()) {
    residues.push_back(residue);
  }
  else {
    AcidPtr acid = residue->getAcidPtr();
    ResiduePtr mut_residue = getResiduePtrByAcidPtm(residue_list, acid, ptm);
    if (mut_residue.get() == nullptr) {
      LOG4CXX_ERROR(logger, "Proteoform:: residue not found");
      throw("Residue not found");
    }
    residues.push_back(mut_residue);
    change_list.push_back(ChangePtr(new Change(0,1, PROTEIN_VARIABLE_CHANGE, 
                                               ptm->getMonoMass(), ptm)));
  }
  for (int i = len + 1; i < residue_seq_ptr_->getLen(); i++) {
    residues.push_back(residue_seq_ptr_->getResiduePtr(i));
  }
  std::string name = name_ + " " + prot_mod_ptr->getName();
  ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(residues));
  for (int i = 0; i < seq_ptr->getLen(); i++) {
    PtmPtr ptm_ptr = seq_ptr->getResiduePtr(i)->getPtmPtr();
    if (!ptm_ptr->isEmpty()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(i, i+1, FIXED_CHANGE, ptm_ptr->getMonoMass(), ptm_ptr));
      change_list.push_back(change_ptr);
    }
  }
  return ProteoformPtr(
      new Proteoform(ori_form_ptr, name, seq_ptr, len, residue_seq_ptr_->getLen()-1, change_list));
}

std::string Proteoform::toString() {
  std::stringstream s;
  s<< "Name: " << name_ << std::endl;
  s<< "Begin pos: " << start_pos_ << std::endl;
  s<< "End pos: " << end_pos_ << std::endl;
  s<< "String: " << residue_seq_ptr_->toString();
  return s.str();

}

ProteoformPtr getOriProteoformPtr(std::string name, ResSeqPtr res_seq_ptr) {
  int start_pos = 0;
  int end_pos = res_seq_ptr->getLen() - 1;
  ChangePtrVec change_list;  
  for (int i = 0; i < res_seq_ptr->getLen(); i++) {
    PtmPtr ptm_ptr = res_seq_ptr->getResiduePtr(i)->getPtmPtr();
    if (!ptm_ptr->isEmpty()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(i, i+1, FIXED_CHANGE, ptm_ptr->getMonoMass(), ptm_ptr));
      change_list.push_back(change_ptr);
    }
  }
  ProteoformPtr ptr = ProteoformPtr(nullptr);
  return ProteoformPtr(new Proteoform(ptr, name, res_seq_ptr, start_pos, end_pos, change_list));
}

ProteoformPtrVec generateProtModProteoform(ProteoformPtrVec &ori_forms, 
                                           ResiduePtrVec &residue_list,
                                           ProtModPtrVec &prot_mods) {
  ProteoformPtrVec new_forms;
  for (unsigned int i = 0; i < ori_forms.size(); i++) {
    for (unsigned int j = 0; j < prot_mods.size(); j++) {
      ProteoformPtr ptr = ori_forms[i]->getProtModProteoform(ori_forms[i], residue_list, prot_mods[j]);
      if (ptr.get() != nullptr) {
        new_forms.push_back(ptr);
      }
    }
  }
  return new_forms;
}

} /* namespace prot */

