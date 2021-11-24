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

#include "common/util/logger.hpp"
#include "common/base/amino_acid_base.hpp"
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
  ProtModPtrVec prot_mod_list = mng_ptr->prsm_para_ptr_->getProtModPtrVec();
  ProtModPtr m_actyl = ProtModBase::getProtModPtr_M_ACETYLATION();
  //LOG_DEBUG("form start " << form_start_pos);
  if (form_start_pos == 0) {
    if (prot_mod_util::containMod(prot_mod_list, m_actyl) && 
        prot_mod_util::allowMod(m_actyl, db_res_seq_ptr->getResidues())) {
      //LOG_DEBUG("M aceylation modifiable!");
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
    bool found = false;
    for (size_t t = 0; t < prot_mod_list.size(); t++) {
      ProtModPtr cur_mod_ptr = prot_mod_list[t];
      if (cur_mod_ptr->getType() == ProtModBase::getType_NME_ACETYLATION() 
          && prot_mod_util::allowMod(cur_mod_ptr, db_res_seq_ptr->getResidues())) {
        found = true;
        prot_mod_ptr = cur_mod_ptr;
        break;
      }
    }
    if (found) {
      //LOG_DEBUG("NME aceylation modifiable!");
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

ProteoformPtr createProteoformPtr(ProteoformPtr base_form_ptr, int match_score, 
                                  std::vector<double> scr_vec, PtmPtr ptm_ptr,
                                  LocalMngPtr mng_ptr) {
  int bgn, end;
  double conf;
  local_util::scrFilter(scr_vec, bgn, end, conf, mng_ptr->threshold_);

  if (bgn == -1) return nullptr;

  //LOG_DEBUG("local start " << bgn << " end " << end);
  LocalAnnoPtr anno = std::make_shared<LocalAnno>(bgn, end, conf, scr_vec,
                                                  match_score, ptm_ptr);
  AminoAcidPtr empty_aa_ptr = AminoAcidBase::getEmptyAminoAcidPtr();
  ModPtr mod_ptr = std::make_shared<Mod>(ResidueBase::getEmptyResiduePtr(),
                                         ResidueBase::getBaseResiduePtr(empty_aa_ptr, ptm_ptr));

  AlterPtr alter = std::make_shared<Alter>(anno->getLeftBpPos(),
                                           anno->getRightBpPos() + 1,
                                           AlterType::VARIABLE, 
                                           ptm_ptr->getMonoMass(),  
                                           mod_ptr);

  alter->setLocalAnno(anno);

  MassShiftPtr mass_shift_ptr = std::make_shared<MassShift>(alter);

  MassShiftPtrVec all_shift_ptr_vec = base_form_ptr->getMassShiftPtrVec(); 
  all_shift_ptr_vec.push_back(mass_shift_ptr);
  std::sort(all_shift_ptr_vec.begin(), all_shift_ptr_vec.end(), MassShift::cmpPosInc);

  ProteoformPtr shift_form_ptr = std::make_shared<Proteoform>(base_form_ptr->getFastaSeqPtr(),
                                                              base_form_ptr->getProtModPtr(),
                                                              base_form_ptr->getStartPos(),
                                                              base_form_ptr->getEndPos(),
                                                              base_form_ptr->getResSeqPtr(),
                                                              all_shift_ptr_vec);


  return shift_form_ptr;
}

ProteoformPtr createProteoformPtr(ProteoformPtr base_form_ptr, 
                                  int match_score, int break_pos,
                                  std::vector<double> scr_vec_1, 
                                  std::vector<double> scr_vec_2,
                                  PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2,
                                  LocalMngPtr mng_ptr) {
  int bgn_1, end_1;
  double conf_1;
  local_util::scrFilter(scr_vec_1, bgn_1, end_1, conf_1, mng_ptr->threshold_);
  if (bgn_1 == -1) return nullptr;

  int bgn_2, end_2;
  double conf_2;
  local_util::scrFilter(scr_vec_2, bgn_2, end_2, conf_2, mng_ptr->threshold_);
  if (bgn_2 == -1) return nullptr;

  if (end_1 > break_pos || bgn_2 <= break_pos) {
    LOG_ERROR("ERROR in localization!");
    exit(EXIT_FAILURE);
  }

  LocalAnnoPtr anno_1 = std::make_shared<LocalAnno>(bgn_1, end_1, conf_1, scr_vec_1,
                                                    match_score, ptm_ptr_1);

  AminoAcidPtr empty_aa_ptr = AminoAcidBase::getEmptyAminoAcidPtr();
  ModPtr mod_ptr_1 = std::make_shared<Mod>(ResidueBase::getEmptyResiduePtr(),
                                           ResidueBase::getBaseResiduePtr(empty_aa_ptr, ptm_ptr_1));

  AlterPtr alter_1 = std::make_shared<Alter>(anno_1->getLeftBpPos(),
                                             break_pos + 1,
                                             AlterType::VARIABLE, 
                                             ptm_ptr_1->getMonoMass(), 
                                             mod_ptr_1);

  alter_1->setLocalAnno(anno_1);

  LocalAnnoPtr anno_2 = std::make_shared<LocalAnno>(bgn_2, end_2, conf_2, scr_vec_2,
                                                    match_score, ptm_ptr_2);

  ModPtr mod_ptr_2 = std::make_shared<Mod>(ResidueBase::getEmptyResiduePtr(),
                                           ResidueBase::getBaseResiduePtr(empty_aa_ptr, ptm_ptr_2));
  AlterPtr alter_2 = std::make_shared<Alter>(break_pos + 1, 
                                             anno_2->getRightBpPos() + 1,
                                             AlterType::VARIABLE, 
                                             ptm_ptr_2->getMonoMass(),
                                             mod_ptr_2);

  alter_2->setLocalAnno(anno_2);


  MassShiftPtrVec all_shift_ptr_vec = base_form_ptr->getMassShiftPtrVec(); 
  MassShiftPtr unexp_shift_ptr_1 = std::make_shared<MassShift>(alter_1);
  all_shift_ptr_vec.push_back(unexp_shift_ptr_1);
  MassShiftPtr unexp_shift_ptr_2 = std::make_shared<MassShift>(alter_2);
  all_shift_ptr_vec.push_back(unexp_shift_ptr_2);

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

