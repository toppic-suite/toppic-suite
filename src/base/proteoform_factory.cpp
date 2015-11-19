#include <sstream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/ptm_base.hpp"
#include "base/residue_base.hpp"
#include "base/trunc_util.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod_base.hpp"
#include "base/prot_mod_util.hpp"
#include "base/proteoform_factory.hpp"

namespace prot {

ProteoformPtr ProteoformFactory::geneDbProteoformPtr(DbResSeqPtr db_res_seq_ptr) {
  int start_pos = 0;
  int end_pos = db_res_seq_ptr->getLen() - 1;
  ChangePtrVec change_list;
  // change list need to be updated
  /*
  for (int i = 0; i < db_res_seq_ptr->getLen(); i++) {
    PtmPtr ptm_ptr = db_res_seq_ptr->getResiduePtr(i)->getPtmPtr();
    if (ptm_ptr == PtmBase::getEmptyPtmPtr()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(i, i+1, ChangeType::FIXED, ptm_ptr->getMonoMass(), ptm_ptr));
      change_list.push_back(change_ptr);
    }
  }
  */
  ProtModPtr none_prot_mod_ptr = ProtModBase::getProtModPtr_NONE();
  return ProteoformPtr(new Proteoform(db_res_seq_ptr, none_prot_mod_ptr,
                                      db_res_seq_ptr, start_pos, end_pos,
                                      change_list));
}

ProteoformPtr ProteoformFactory::geneProtModProteoform(ProteoformPtr db_form_ptr,
                                                       ProtModPtr prot_mod_ptr) {
  // check if the proteoform can be truncated
  DbResSeqPtr db_res_seq_ptr = db_form_ptr->getDbResSeqPtr();
  bool valid_mod = ProtModUtil::allowMod(prot_mod_ptr, db_res_seq_ptr->getResidues());
  if (!valid_mod) {
    //LOG_DEBUG("NO valid mod");
    return ProteoformPtr(nullptr);
  }

  TruncPtr trunc_ptr = prot_mod_ptr->getTruncPtr();
  int start = trunc_ptr->getTruncLen();
  /* last bp index */
  int end = db_form_ptr->getLen();
  // copy input changes
  ChangePtrVec ori_change_ptrs = db_form_ptr->getChangePtrVec();
  ChangePtrVec change_ptrs;
  for (size_t i = 0; i < ori_change_ptrs.size(); i++) {
    if (ori_change_ptrs[i]->getLeftBpPos() >= start
        && ori_change_ptrs[i]->getRightBpPos() <= end + 1) {
      ChangePtr change_ptr = Change::geneChangePtr(ori_change_ptrs[i], start);
      change_ptrs.push_back(change_ptr);
    }
  }

  ResiduePtrVec ori_vec = db_res_seq_ptr->getResidues();
  ResiduePtrVec new_vec(ori_vec.begin()+start, ori_vec.begin()+end);
  // apply mod
  ModPtr mod_ptr = prot_mod_ptr->getModPtr();
  if (!ModBase::isNoneModPtr(mod_ptr)) {
    int mod_pos = prot_mod_ptr->getModPos();
    int new_pos = mod_pos - start;
    new_vec[new_pos] = prot_mod_ptr->getModPtr()->getModResiduePtr();
    change_ptrs.push_back(ChangePtr(new Change(new_pos, new_pos+1, ChangeType::PROTEIN_VARIABLE,
                                               mod_ptr->getShift(), mod_ptr)));
  }
  ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(new_vec));

  //LOG_DEBUG("mod protein sequence name " << db_res_seq_ptr->getName()
  //<< " len " << db_res_seq_ptr->getLen());
  return ProteoformPtr(
      new Proteoform(db_res_seq_ptr, prot_mod_ptr, seq_ptr, start,
                     db_res_seq_ptr->getLen()-1, change_ptrs));
}

ProteoformPtr ProteoformFactory::geneSubProteoform(ProteoformPtr proteoform_ptr,
                                                   int local_start, int local_end) {
    ResiduePtrVec residues;
    ResSeqPtr res_seq_ptr = proteoform_ptr->getResSeqPtr();
    for (int i = local_start; i <= local_end; i++) {
        residues.push_back(res_seq_ptr->getResiduePtr(i));
    }
    ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(residues));
    ChangePtrVec change_list;
    ChangePtrVec ori_change_list = proteoform_ptr->getChangePtrVec();
    for (size_t i = 0; i < ori_change_list.size(); i++) {
        if (ori_change_list[i]->getLeftBpPos() >= local_start
                && ori_change_list[i]->getRightBpPos() <= local_end + 1) {
            ChangePtr change_ptr = Change::geneChangePtr(ori_change_list[i], local_start);
            change_list.push_back(change_ptr);
        }
    }
    DbResSeqPtr db_res_seq_ptr = proteoform_ptr->getDbResSeqPtr();
    ProtModPtr prot_mod_ptr = proteoform_ptr->getProtModPtr();
    return ProteoformPtr(
               new Proteoform(db_res_seq_ptr, prot_mod_ptr, seq_ptr,
                              local_start + proteoform_ptr->getStartPos(),
                              local_end + proteoform_ptr->getStartPos(), change_list));
}

ProteoformPtrVec ProteoformFactory::geneProtModProteoform(ProteoformPtr proteo_ptr,
                                                          const ProtModPtrVec &prot_mods) {
  ProteoformPtrVec new_forms;
  for (size_t j = 0; j < prot_mods.size(); j++) {
    ProteoformPtr ptr = geneProtModProteoform(proteo_ptr, prot_mods[j]);
    if (ptr.get() != nullptr) {
      new_forms.push_back(ptr);
    }
  }
  return new_forms;
}

ProteoformPtrVec ProteoformFactory::geneProtModProteoform(const ProteoformPtrVec &ori_forms,
                                                          const ProtModPtrVec &prot_mods) {
  ProteoformPtrVec new_forms;
  for (size_t i = 0; i < ori_forms.size(); i++) {
    for (size_t j = 0; j < prot_mods.size(); j++) {
      ProteoformPtr ptr = geneProtModProteoform(ori_forms[i], prot_mods[j]);
      if (ptr.get() != nullptr) {
        new_forms.push_back(ptr);
      }
    }
  }
  return new_forms;
}

ProteoformPtrVec2D ProteoformFactory::gene2DProtModProteoform(const ProteoformPtrVec &ori_forms,
                                                              const ProtModPtrVec &prot_mods) {
  ProteoformPtrVec2D new_forms;
  for (size_t i = 0; i < ori_forms.size(); i++) {
    ProteoformPtrVec mod_forms;
    for (size_t j = 0; j < prot_mods.size(); j++) {
      ProteoformPtr ptr = geneProtModProteoform(ori_forms[i], prot_mods[j]);
      if (ptr.get() != nullptr) {
        mod_forms.push_back(ptr);
      }
    }
    new_forms.push_back(mod_forms);
  }
  return new_forms;
}

} /* namespace prot */

