#include <sstream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/proteoform.hpp"

namespace prot {

Proteoform::Proteoform(ProteoformPtr raw_form_ptr, ProtModPtr prot_mod_ptr, std::string name, 
                       ResSeqPtr res_seq_ptr, int start_pos, int end_pos, ChangePtrVec change_list) {
  raw_form_ptr_ = raw_form_ptr;
  prot_mod_ptr_ = prot_mod_ptr;
  name_ = name;
  residue_seq_ptr_ = res_seq_ptr;
  start_pos_ = start_pos;
  end_pos_ = end_pos;
  LOG_TRACE( "start bp spec generation");
  bp_spec_ptr_ = BpSpecPtr(new BpSpec(res_seq_ptr));
  change_list_ = change_list;
  std::sort(change_list.begin(), change_list.end(), compareChangeUp);
}

SegmentPtrVec Proteoform::getSegmentPtrVec() {
  ChangePtrVec unexpected_changes;
  double mass_shift_sum = 0;
  for (unsigned int i = 0; i < change_list_.size(); i++) {
    if (change_list_[i]->getChangeType() == UNEXPECTED_CHANGE) {
      unexpected_changes.push_back(change_list_[i]);
      mass_shift_sum += change_list_[i]->getMassShift();
    }
  }
  SegmentPtrVec segments;
  double n_shift = 0;
  double c_shift = mass_shift_sum;
  int left = 0;
  for (unsigned int i = 0; i < unexpected_changes.size(); i++) {
    int right = unexpected_changes[i]->getLeftBpPos();
    SegmentPtr segment_ptr = SegmentPtr(
        new Segment(left, right, n_shift, c_shift)); 
    segments.push_back(segment_ptr);
    left = unexpected_changes[i]->getRightBpPos();
    n_shift = n_shift + unexpected_changes[i]->getMassShift();
    c_shift = c_shift - unexpected_changes[i]->getMassShift();
  }
  int right = residue_seq_ptr_->getLen();
  SegmentPtr segment_ptr = SegmentPtr(
      new Segment(left, right, n_shift, c_shift)); 
  segments.push_back(segment_ptr);
  return segments;
}

bool isValidTrunc(ProteoformPtr raw_form_ptr, ProtModPtr prot_mod_ptr) {
  TruncPtr trunc_ptr = prot_mod_ptr->getTruncPtr();
  int trunc_len = trunc_ptr->getTruncLen();
  AcidPtrVec trunc_acids = trunc_ptr->getAcidPtrVec();
  ResSeqPtr res_seq_ptr = raw_form_ptr->getResSeqPtr();  
  //check if trunc acids match N-terminal acids of the protein 
  if (trunc_len >= res_seq_ptr->getLen()) {
    return false; ;
  }
  bool match = true;
  for (int i = 0; i < trunc_len; i++) {
    AcidPtr acid = res_seq_ptr->getResiduePtr(i)->getAcidPtr();
    if (acid.get() != trunc_acids[i].get()) {
      match = false;
      break;
    }
  }
  return match;
}

ProteoformPtr getProtModProteoform(ProteoformPtr raw_form_ptr,
                                   ResiduePtrVec &residue_list, 
                                   ProtModPtr prot_mod_ptr) {
  bool valid_trunc = isValidTrunc(raw_form_ptr, prot_mod_ptr);
  if (!valid_trunc) {
    return ProteoformPtr(nullptr);
  }
  // first residue might be acetylated 
  ResSeqPtr res_seq_ptr = raw_form_ptr->getResSeqPtr();  
  int start = prot_mod_ptr->getTruncPtr()->getTruncLen();
  ResiduePtrVec residues;
  ResiduePtr residue = res_seq_ptr->getResiduePtr(start);
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
      LOG_ERROR( "Proteoform:: residue not found");
      throw("Residue not found");
    }
    residues.push_back(mut_residue);
    change_list.push_back(ChangePtr(new Change(0,1, PROTEIN_VARIABLE_CHANGE, 
                                               ptm->getMonoMass(), ptm)));
  }
  // add all other residues
  for (int i = start + 1; i < res_seq_ptr->getLen(); i++) {
    residues.push_back(res_seq_ptr->getResiduePtr(i));
  }
  ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(residues));
  for (int i = 0; i < seq_ptr->getLen(); i++) {
    PtmPtr ptm_ptr = seq_ptr->getResiduePtr(i)->getPtmPtr();
    if (!ptm_ptr->isEmpty()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(i, i+1, FIXED_CHANGE, ptm_ptr->getMonoMass(), ptm_ptr));
      change_list.push_back(change_ptr);
    }
  }
  std::string name = raw_form_ptr->getName() + " " + prot_mod_ptr->getName();
  return ProteoformPtr(
      new Proteoform(raw_form_ptr, prot_mod_ptr, name, seq_ptr, start, 
                     res_seq_ptr->getLen()-1, change_list));
}

ProteoformPtr getSubProteoform(ProteoformPtr proteoform_ptr, int start, int end) {
  ResiduePtrVec residues;
  ResSeqPtr res_seq_ptr = proteoform_ptr->getResSeqPtr();  
  for (int i = start; i <= end; i++) {
    residues.push_back(res_seq_ptr->getResiduePtr(i));
  }
  ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(residues));
  ChangePtrVec change_list;
  ChangePtrVec ori_change_list = proteoform_ptr->getChangePtrVec();
  for (unsigned int i = 0; i < ori_change_list.size(); i++) {
    if (ori_change_list[i]->getLeftBpPos() >= start 
        && ori_change_list[i]->getRightBpPos() <= end + 1) {
      ChangePtr change_ptr = ChangePtr(new Change(*ori_change_list[i], start));
      change_list.push_back(change_ptr);
    }
  }
  ProteoformPtr raw_form_ptr = proteoform_ptr->getRawFormPtr();
  ProtModPtr prot_mod_ptr = proteoform_ptr->getProtModPtr();
  std::string name = proteoform_ptr->getName();

  return ProteoformPtr(
      new Proteoform(raw_form_ptr, prot_mod_ptr, name, seq_ptr, start + proteoform_ptr->getStartPos(), 
                     end + proteoform_ptr->getEndPos(), change_list));
}

std::string Proteoform::toString() {
  std::stringstream s;
  s<< "Name: " << name_ << std::endl;
  s<< "Begin pos: " << start_pos_ << std::endl;
  s<< "End pos: " << end_pos_ << std::endl;
  s<< "String: " << residue_seq_ptr_->toString();
  return s.str();
}

ProteoformPtr getRawProteoformPtr(std::string name, ResSeqPtr res_seq_ptr) {
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
  ProteoformPtr proteoform_ptr = ProteoformPtr(nullptr);
  ProtModPtr prot_mod_ptr = ProtModPtr(nullptr);
  return ProteoformPtr(new Proteoform(proteoform_ptr, prot_mod_ptr, name, 
                                      res_seq_ptr, start_pos, end_pos, change_list));
}

ProteoformPtrVec generateProtModProteoform(ProteoformPtrVec &ori_forms, 
                                           ResiduePtrVec &residue_list,
                                           ProtModPtrVec &prot_mods) {
  ProteoformPtrVec new_forms;
  for (unsigned int i = 0; i < ori_forms.size(); i++) {
    for (unsigned int j = 0; j < prot_mods.size(); j++) {
      ProteoformPtr ptr = getProtModProteoform(ori_forms[i], residue_list, prot_mods[j]);
      if (ptr.get() != nullptr) {
        new_forms.push_back(ptr);
      }
    }
  }
  return new_forms;
}

ProteoformPtrVec getProtModNoneProteoform(ProteoformPtrVec &all_forms) {
  ProteoformPtrVec new_forms;
  for (unsigned int i = 0; i < all_forms.size(); i++) {
    ProtModPtr prot_mod_ptr = all_forms[i]->getProtModPtr();
    if (prot_mod_ptr.get() != nullptr && prot_mod_ptr->isNone()) {
      new_forms.push_back(all_forms[i]);
    }
  }
  return new_forms;
}


} /* namespace prot */

