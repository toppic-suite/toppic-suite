//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <algorithm>

#include "common/base/residue_base.hpp"
#include "common/base/prot_mod_util.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/residue_util.hpp"

#include "stat/local/local_util.hpp"
#include "stat/local/local_proteoform.hpp"

namespace toppic {

namespace local_proteoform {

void getNtermTruncRange(ProteoformPtr proteoform, LocalMngPtr mng_ptr, int & min, int & max) {

  int ori_start = proteoform->getStartPos();
  int ori_end = proteoform->getEndPos();

  min = ori_start - mng_ptr->n_term_range_;
  if (min < 0) {
    min = 0;
  }

  max = ori_start + mng_ptr->n_term_range_;
  if (max > ori_end) {
    max = ori_end;
  }
}

void getCtermTruncRange(ProteoformPtr proteoform, LocalMngPtr mng_ptr, int & min, int & max) {

  int ori_start = proteoform->getStartPos();
  int ori_end = proteoform->getEndPos();

  min = ori_end - mng_ptr->c_term_range_;
  if (min < ori_start) {
    min = ori_start;
  }

  int raw_seq_len = static_cast<int>(proteoform->getFastaSeqPtr()->getRawSeq().length());
  max = ori_end + mng_ptr->c_term_range_;
  if (max > raw_seq_len - 1) {
    max = raw_seq_len - 1;
  }
}

ProteoformPtrVec getCandidateForm(FastaSeqPtr seq_ptr, int ori_start_pos, 
                                  int form_start_pos, int form_end_pos, 
                                  MassShiftPtrVec & exp_shift_ptr_vec, 
                                  LocalMngPtr mng_ptr) {
  ProteoformPtrVec result_vec;
  // update mass shifts
  MassShiftPtrVec valid_shift_ptr_vec;
  for (size_t i = 0; i < exp_shift_ptr_vec.size(); i++) {
    MassShiftPtr shift_ptr = exp_shift_ptr_vec[i];
    int shift_start_pos = shift_ptr->getLeftBpPos() + ori_start_pos;
    int shift_end_pos = shift_ptr->getRightBpPos() + ori_start_pos - 1;
    if (shift_start_pos >= form_start_pos && shift_end_pos <= form_end_pos) {
      int left_bp_pos = shift_ptr->getLeftBpPos() + ori_start_pos - form_start_pos;
      int right_bp_pos = shift_ptr->getRightBpPos() + ori_start_pos - form_start_pos;
      AlterPtr alter_ptr = std::make_shared<Alter>(left_bp_pos, right_bp_pos,  
                                                   shift_ptr->getTypePtr(), 
                                                   shift_ptr->getMassShift(), 
                                                   shift_ptr->getAlterPtr(0)->getModPtr());
      MassShiftPtr new_shift_ptr = std::make_shared<MassShift>(alter_ptr);
      valid_shift_ptr_vec.push_back(new_shift_ptr);
    }
  }

  // obtain res_seq
  std::string whole_seq = seq_ptr->getRawSeq();
  std::string sub_seq = whole_seq.substr(form_start_pos, form_end_pos - form_start_pos + 1); 
  ModPtrVec fixed_mod_ptr_vec = mng_ptr->prsm_para_ptr_->getFixModPtrVec();
  ResSeqPtr db_res_seq_ptr = std::make_shared<ResidueSeq>(residue_util::convertStrToResiduePtrVec(whole_seq, fixed_mod_ptr_vec));
  ResSeqPtr sub_res_seq_ptr = std::make_shared<ResidueSeq>(residue_util::convertStrToResiduePtrVec(sub_seq, fixed_mod_ptr_vec));

  // obtain protein mod
  ProtModPtr prot_mod_ptr = ProtModBase::getProtModPtr_NONE(); 
  if (form_start_pos == 1) {
    ProtModPtr nme = ProtModBase::getProtModPtr_NME();
    if (prot_mod_util::containMod(mng_ptr->prsm_para_ptr_->getProtModPtrVec(), nme) && 
        prot_mod_util::allowMod(nme, db_res_seq_ptr->getResidues())) {
      prot_mod_ptr = nme;
    }
  }

  ProteoformPtr form_ptr
    = std::make_shared<Proteoform>(seq_ptr, 
                                   prot_mod_ptr, 
                                   form_start_pos, 
                                   form_end_pos, 
                                   sub_res_seq_ptr,
                                   valid_shift_ptr_vec);
  result_vec.push_back(form_ptr);

  // check possible N-terminal acetylation form
  if (form_start_pos == 0) {
    ProtModPtr m_actyl = ProtModBase::getProtModPtr_M_ACETYLATION();
    if (prot_mod_util::containMod(mng_ptr->prsm_para_ptr_->getProtModPtrVec(), m_actyl) && 
        prot_mod_util::allowMod(m_actyl, db_res_seq_ptr->getResidues())) {
      prot_mod_ptr = m_actyl;
      // apply mod
      ModPtr mod_ptr = prot_mod_ptr->getModPtr();
      ResiduePtrVec res_vec = sub_res_seq_ptr->getResidues();
      int pos = 0;
      res_vec[pos] = mod_ptr->getModResiduePtr();
      ResSeqPtr new_seq_ptr = std::make_shared<ResidueSeq>(res_vec);

      AlterPtr alter_ptr = std::make_shared<Alter>(pos, pos + 1, 
                                                   AlterType::PROTEIN_VARIABLE,
                                                   mod_ptr->getShift(), mod_ptr);
      MassShiftPtr shift_ptr = std::make_shared<MassShift>(alter_ptr); 
      valid_shift_ptr_vec.push_back(shift_ptr);

      form_ptr = std::make_shared<Proteoform>(seq_ptr, 
                                              prot_mod_ptr, 
                                              form_start_pos, 
                                              form_end_pos, 
                                              new_seq_ptr,
                                              valid_shift_ptr_vec);
      result_vec.push_back(form_ptr);
    }
  }
  else if (form_start_pos == 1) {
    ProtModPtr nme_actyl = ProtModBase::getProtModPtr_NME_ACETYLATION();
    if (prot_mod_util::containMod(mng_ptr->prsm_para_ptr_->getProtModPtrVec(), nme_actyl) && 
        prot_mod_util::allowMod(nme_actyl, db_res_seq_ptr->getResidues())) {
      prot_mod_ptr = nme_actyl;
      // apply mod
      ModPtr mod_ptr = prot_mod_ptr->getModPtr();
      ResiduePtrVec res_vec = sub_res_seq_ptr->getResidues();
      int pos = 0;
      res_vec[pos] = mod_ptr->getModResiduePtr();
      ResSeqPtr new_seq_ptr = std::make_shared<ResidueSeq>(res_vec);

      AlterPtr alter_ptr = std::make_shared<Alter>(pos, pos + 1, 
                                                   AlterType::PROTEIN_VARIABLE,
                                                   mod_ptr->getShift(), mod_ptr);
      MassShiftPtr shift_ptr = std::make_shared<MassShift>(alter_ptr); 
      valid_shift_ptr_vec.push_back(shift_ptr);

      form_ptr = std::make_shared<Proteoform>(seq_ptr, 
                                              prot_mod_ptr, 
                                              form_start_pos, 
                                              form_end_pos, 
                                              new_seq_ptr,
                                              valid_shift_ptr_vec);
      result_vec.push_back(form_ptr);
    }
  }
  return result_vec;
}


ProteoformPtrVec getAllCandidateForms(ProteoformPtr ori_form_ptr, 
                                      LocalMngPtr mng_ptr) {
  ProteoformPtrVec result_form_vec; 

  MassShiftPtrVec ori_shift_ptr_vec = ori_form_ptr->getMassShiftPtrVec();
  // Filter out unexpected mass shifts and protein N-terminal variable PTMs
  MassShiftPtrVec tmp_vec = local_util::massShiftFilter(ori_shift_ptr_vec,
                                                        AlterType::UNEXPECTED);
  MassShiftPtrVec exp_shift_ptr_vec = local_util::massShiftFilter(tmp_vec, 
                                                                  AlterType::PROTEIN_VARIABLE);
  FastaSeqPtr fasta_seq_ptr = ori_form_ptr->getFastaSeqPtr();
  int ori_start_pos = ori_form_ptr->getStartPos();
  int ori_end_pos = ori_form_ptr->getEndPos();
  // 1. Add the orignal form
  ProteoformPtrVec base_form_ptr_vec = getCandidateForm(fasta_seq_ptr, 
                                                        ori_start_pos,
                                                        ori_start_pos,
                                                        ori_end_pos,
                                                        exp_shift_ptr_vec,
                                                        mng_ptr);
  result_form_vec.insert(std::end(result_form_vec), std::begin(base_form_ptr_vec), std::end(base_form_ptr_vec));

  // 2. N-term trancations
  int start_pos_min, start_pos_max; 
  getNtermTruncRange(ori_form_ptr, mng_ptr, start_pos_min, start_pos_max);
  for (int i = start_pos_min; i <= start_pos_max; i++) {
    if (i == ori_start_pos) {
      continue;
    }
    ProteoformPtrVec form_ptr_vec = getCandidateForm(fasta_seq_ptr, 
                                                     ori_start_pos,
                                                     i,
                                                     ori_end_pos,
                                                     exp_shift_ptr_vec,
                                                     mng_ptr);
    result_form_vec.insert(std::end(result_form_vec), std::begin(form_ptr_vec), std::end(form_ptr_vec));
  }

  // 3. C-term truncations
  int end_pos_min, end_pos_max; 
  getCtermTruncRange(ori_form_ptr, mng_ptr, end_pos_min, end_pos_max);
  for (int i = end_pos_min; i <= end_pos_max; i++) {
    if (i == ori_end_pos) {
      continue;
    }
    ProteoformPtrVec form_ptr_vec = getCandidateForm(fasta_seq_ptr, 
                                                     ori_start_pos,
                                                     ori_start_pos,
                                                     i,
                                                     exp_shift_ptr_vec, 
                                                     mng_ptr);
    result_form_vec.insert(std::end(result_form_vec), std::begin(form_ptr_vec), std::end(form_ptr_vec));
  }

  return result_form_vec;
}

ProteoformPtr createProteoformPtr(ProteoformPtr base_form_ptr, 
                                  double shift_mass, int match_score, 
                                  std::vector<double> scr_vec, PtmPtr ptm_ptr,
                                  LocalMngPtr mng_ptr) {
  int bgn, end;
  double conf;
  local_util::scrFilter(scr_vec, bgn, end, conf, mng_ptr->threshold_);

  if (bgn == -1) return nullptr;

  LocalAnnoPtr anno = std::make_shared<LocalAnno>(bgn, end, conf, scr_vec,
                                                  match_score, ptm_ptr);

  AlterPtr alter = std::make_shared<Alter>(anno->getLeftBpPos(),
                                           anno->getRightBpPos() + 1,
                                           AlterType::UNEXPECTED, 
                                           shift_mass,
                                           std::make_shared<Mod>(ResidueBase::getEmptyResiduePtr(),
                                                                 ResidueBase::getEmptyResiduePtr()));

  alter->setLocalAnno(anno);

  MassShiftPtr unexp_shift_ptr = std::make_shared<MassShift>(alter);

  MassShiftPtrVec all_shift_ptr_vec = base_form_ptr->getMassShiftPtrVec(); 
  all_shift_ptr_vec.push_back(unexp_shift_ptr);
  std::sort(all_shift_ptr_vec.begin(), all_shift_ptr_vec.end(), MassShift::cmpPosInc);

  ProteoformPtr shift_form_ptr = std::make_shared<Proteoform>(base_form_ptr->getFastaSeqPtr(),
                                                              base_form_ptr->getProtModPtr(),
                                                              base_form_ptr->getStartPos(),
                                                              base_form_ptr->getEndPos(),
                                                              base_form_ptr->getResSeqPtr(),
                                                              all_shift_ptr_vec);


  return shift_form_ptr;
}

}

}

