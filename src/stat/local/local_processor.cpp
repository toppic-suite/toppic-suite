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
#include <cmath>
#include <numeric>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/base/ptm.hpp"
#include "common/base/ptm_util.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/residue_util.hpp"
#include "common/base/prot_mod.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/prot_mod_util.hpp"

#include "seq/local_anno.hpp"
#include "seq/mass_shift.hpp"
#include "seq/proteoform_factory.hpp"

#include "ms/spec/msalign_util.hpp"
#include "ms/spec/spectrum_set.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"

#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"

#include "stat/local/local_util.hpp"
#include "stat/local/local_processor.hpp"

namespace toppic {

LocalProcessor::LocalProcessor(LocalMngPtr mng_ptr):
  mng_ptr_(mng_ptr) {
    init();
  }

void LocalProcessor::init() {
  ptm_ptr_vec_ = ptm_util::readPtmTxt(mng_ptr_->residueModFileName_);

  for (size_t i = 0; i < ptm_ptr_vec_.size(); i++) {
    for (size_t j = 0; j < ptm_ptr_vec_.size(); j++) {
      ptm_pair_vec_.push_back(std::make_pair(ptm_ptr_vec_[i], ptm_ptr_vec_[j]));
    }
  }

  std::vector<ModPtrVec> mod_ptr_vec2d = mod_util::readModTxt(mng_ptr_->residueModFileName_);
  mod_list_N_ = mod_ptr_vec2d[0];
  mod_list_C_ = mod_ptr_vec2d[1];
  mod_list_any_ = mod_ptr_vec2d[2];
}

void LocalProcessor::process() {
  std::string spec_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();

  std::string input_file_name = file_util::basename(spec_file_name) + "." + mng_ptr_->input_file_ext_;

  std::string output_file_name = file_util::basename(spec_file_name) + "." + mng_ptr_->output_file_ext_;

  PrsmXmlWriterPtr prsm_writer = std::make_shared<PrsmXmlWriter>(output_file_name);

  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();

  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);

  PrsmReaderPtr prsm_reader = std::make_shared<PrsmReader>(input_file_name);

  ModPtrVec fix_mod_list = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();

  PrsmPtr prsm_ptr = prsm_reader->readOnePrsm(seq_reader, fix_mod_list);

  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();

  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();

  SimpleMsAlignReaderPtr reader_ptr = std::make_shared<SimpleMsAlignReader>(spec_file_name, 
                                                                            group_spec_num,
                                                                            sp_para_ptr->getActivationPtr());

  SpectrumSetPtr spec_set_ptr;

  int spectrum_num = msalign_util::getSpNum(mng_ptr_->prsm_para_ptr_->getSpectrumFileName());

  int cnt = 0;

  while ((spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(reader_ptr, sp_para_ptr))!= nullptr) {
    cnt += group_spec_num;
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec
            = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);

        if (prsm_ptr->getProteoformPtr()->getMassShiftNum(AlterType::UNEXPECTED) > 0) {
          prsm_ptr = processOnePrsm(prsm_ptr);
        }

        prsm_writer->write(prsm_ptr);
        prsm_ptr = prsm_reader->readOnePrsm(seq_reader, fix_mod_list);
      }
    }
    std::cout << std::flush << "PTM characterization is processing " << cnt
        << " of " << spectrum_num << " spectra.\r";
  }
  prsm_reader->close();
  prsm_writer->close();
  std::cout << std::endl;
}

PrsmPtr LocalProcessor::processOnePrsm(PrsmPtr prsm) {
  int mass_shift_num = prsm->getProteoformPtr()->getMassShiftNum(AlterType::UNEXPECTED);
  double err_tole = mng_ptr_->peak_tole_ptr_->compStrictErrorTole(prsm->getAdjustedPrecMass());
  if (mass_shift_num == 1) {
    double mass = prsm->getProteoformPtr()->getMassShiftPtrVec(AlterType::UNEXPECTED)[0]->getMassShift();
    // if the mass shift is between [-1, 1], we don't characterize it
    if (std::abs(mass) <= 1 + err_tole) return prsm;
    return processOnePtm(prsm);
  } else if (mass_shift_num == 2) {
    //TO DO We will add the function of processTwoPtm later.
    //return processTwoPtm(prsm);
  }
  return prsm;
}

PrsmPtr LocalProcessor::processOnePtm(PrsmPtr prsm) {
  int ori_num_match_ion = local_util::compMatchFragNum(prsm->getProteoformPtr(),
                                                       prsm->getRefineMsPtrVec(),
                                                       mng_ptr_->min_mass_);

  // we will get a nullptr if the mass shift can't be explained by a known variable ptm
  ProteoformPtr one_known_proteoform = processOneKnownPtm(prsm);

  if (one_known_proteoform != nullptr) {
    int new_num_match_ion = local_util::compMatchFragNum(one_known_proteoform,
                                                         prsm->getRefineMsPtrVec(),
                                                         mng_ptr_->min_mass_);
    if (new_num_match_ion > ori_num_match_ion - mng_ptr_->DESC_MATCH_LIMIT_
        && new_num_match_ion > ori_num_match_ion * mng_ptr_->desc_ratio_) {
      one_known_proteoform->setProteoClusterId(prsm->getProteoformPtr()->getProteoClusterId());
      one_known_proteoform->setProtId(prsm->getProteoformPtr()->getProtId());
      prsm->setProteoformPtr(one_known_proteoform, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      return prsm;
    }
  }

  /* TO DO we will add two known prsm later
  ProteoformPtr two_known_prsm = processTwoKnownPtm(prsm);

  if (two_known_prsm != nullptr) {
    double new_num_match_ion = local_util::compMatchFragNum(two_known_prsm,
                                                            prsm->getRefineMsPtrVec(),
                                                            mng_ptr_->min_mass_);
    if (new_num_match_ion > ori_num_match_ion - mng_ptr_->DESC_MATCH_LIMIT_
        && new_num_match_ion > ori_num_match_ion * mng_ptr_->desc_ratio_) {
      prsm->setProteoformPtr(two_known_prsm, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      return prsm;
    }
  }
  */

  return prsm;
}

void LocalProcessor::getNtermTruncRange(ProteoformPtr proteoform, int & min, int & max) {

  int ori_start = proteoform->getStartPos();
  int ori_end = proteoform->getEndPos();

  min = ori_start - mng_ptr_->n_term_range_;
  if (min < 0) {
    min = 0;
  }

  max = ori_start + mng_ptr_->n_term_range_;
  if (max > ori_end) {
    max = ori_end;
  }
}

void LocalProcessor::getCtermTruncRange(ProteoformPtr proteoform, int & min, int & max) {

  int ori_start = proteoform->getStartPos();
  int ori_end = proteoform->getEndPos();

  min = ori_end - mng_ptr_->c_term_range_;
  if (min < ori_start) {
    min = ori_start;
  }

  int raw_seq_len = static_cast<int>(proteoform->getFastaSeqPtr()->getRawSeq().length());
  max = ori_end + mng_ptr_->c_term_range_;
  if (max > raw_seq_len - 1) {
    max = raw_seq_len - 1;
  }
}

ProteoformPtrVec LocalProcessor::createCandidateForm(FastaSeqPtr seq_ptr, int ori_start_pos, 
                                                     int form_start_pos, int form_end_pos, 
                                                     MassShiftPtrVec & exp_shift_ptr_vec) {
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
      MassShiftPtr new_shift_ptr = std::make_shared<MassShift>(left_bp_pos, right_bp_pos, shift_ptr->getMassShift());
      valid_shift_ptr_vec.push_back(new_shift_ptr);
    }
  }

  // obtain res_seq
  std::string whole_seq = seq_ptr->getRawSeq();
  std::string sub_seq = whole_seq.substr(form_start_pos, form_end_pos - form_start_pos + 1); 
  ModPtrVec fixed_mod_ptr_vec = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();
  ResSeqPtr db_res_seq_ptr = std::make_shared<ResidueSeq>(residue_util::convertStrToResiduePtrVec(whole_seq, fixed_mod_ptr_vec));
  ResSeqPtr sub_res_seq_ptr = std::make_shared<ResidueSeq>(residue_util::convertStrToResiduePtrVec(sub_seq, fixed_mod_ptr_vec));

  // obtain protein mod
  ProtModPtr prot_mod_ptr = ProtModBase::getProtModPtr_NONE(); 
  if (form_start_pos == 1) {
    ProtModPtr nme = ProtModBase::getProtModPtr_NME();
    if (prot_mod_util::containMod(mng_ptr_->prsm_para_ptr_->getProtModPtrVec(), nme) && 
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
    if (prot_mod_util::containMod(mng_ptr_->prsm_para_ptr_->getProtModPtrVec(), m_actyl) && 
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
    if (prot_mod_util::containMod(mng_ptr_->prsm_para_ptr_->getProtModPtrVec(), nme_actyl) && 
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


ProteoformPtrVec LocalProcessor::getOneKnownPtmCandidateForms(ProteoformPtr ori_form_ptr) { 
  ProteoformPtrVec result_form_vec; 

  MassShiftPtrVec ori_shift_ptr_vec = ori_form_ptr->getMassShiftPtrVec();
  // Filter out unexpected mass shifts and protein N-terminal variable PTMs
  MassShiftPtrVec tmp_vec = local_util::massShiftFilter(ori_shift_ptr_vec,
                                                        AlterType::UNEXPECTED);
  MassShiftPtrVec exp_shift_ptr_vec = local_util::massShiftFilter(tmp_vec, 
                                                                  AlterType::UNEXPECTED);
  FastaSeqPtr fasta_seq_ptr = ori_form_ptr->getFastaSeqPtr();
  int ori_start_pos = ori_form_ptr->getStartPos();
  int ori_end_pos = ori_form_ptr->getEndPos();
  // 1. Add the orignal form
  ProteoformPtrVec base_form_ptr_vec = createCandidateForm(fasta_seq_ptr, 
                                                           ori_start_pos,
                                                           ori_start_pos,
                                                           ori_end_pos,
                                                           exp_shift_ptr_vec);
  result_form_vec.insert(std::end(result_form_vec), std::begin(base_form_ptr_vec), std::end(base_form_ptr_vec));

  // 2. N-term trancations
  int start_pos_min, start_pos_max; 
  getNtermTruncRange(ori_form_ptr, start_pos_min, start_pos_max);
  for (int i = start_pos_min; i <= start_pos_max; i++) {
    if (i == ori_start_pos) {
      continue;
    }
    ProteoformPtrVec form_ptr_vec = createCandidateForm(fasta_seq_ptr, 
                                                        ori_start_pos,
                                                        i,
                                                        ori_end_pos,
                                                        exp_shift_ptr_vec);
    result_form_vec.insert(std::end(result_form_vec), std::begin(form_ptr_vec), std::end(form_ptr_vec));
  }

  // 3. C-term truncations
  int end_pos_min, end_pos_max; 
  getCtermTruncRange(ori_form_ptr, end_pos_min, end_pos_max);
  for (int i = end_pos_min; i <= end_pos_max; i++) {
    if (i == ori_end_pos) {
      continue;
    }
    ProteoformPtrVec form_ptr_vec = createCandidateForm(fasta_seq_ptr, 
                                                        ori_start_pos,
                                                        ori_start_pos,
                                                        i,
                                                        exp_shift_ptr_vec);
    result_form_vec.insert(std::end(result_form_vec), std::begin(form_ptr_vec), std::end(form_ptr_vec));
  }

  return result_form_vec;
}

// we will get a nullptr if the mass shift can't be explained by a variable ptm
ProteoformPtr LocalProcessor::processOneKnownPtm(PrsmPtr prsm_ptr) {

  //get canidate forms
  ProteoformPtr ori_form_ptr = prsm_ptr->getProteoformPtr();
  ProteoformPtrVec cand_form_vec = getOneKnownPtmCandidateForms(ori_form_ptr);

  ExtendMsPtrVec extend_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  double adjust_prec_mass = prsm_ptr->getAdjustedPrecMass();
  double err_tole = mng_ptr_->peak_tole_ptr_->compStrictErrorTole(adjust_prec_mass);
  
  int best_match_score = -1;
  ProteoformPtr best_form_ptr = nullptr;

  for (size_t i = 0; i < cand_form_vec.size(); i++) {
    ProteoformPtr cand_form_ptr = cand_form_vec[i];
    LocalResultPtr local_result_ptr = onePtmLocalize(cand_form_ptr, extend_ms_ptr_vec, 
                                                     adjust_prec_mass, err_tole); 
    if (local_result_ptr != nullptr && local_result_ptr->match_score_ > best_match_score) {
      best_match_score = local_result_ptr->match_score_;
      best_form_ptr = local_result_ptr->form_ptr_;
    }
  }
  return best_form_ptr;
}

LocalResultPtr LocalProcessor::onePtmLocalize(ProteoformPtr form_ptr, const ExtendMsPtrVec & extend_ms_ptr_vec, 
                                              double prec_mass, double err_tole) {

  // Check if the unexpected mass shift can be explained by a PTM
  // sum of all expected shift values, in most cases, there is only one unexpected shift
  double unexp_shift_mass = prec_mass - form_ptr->getMass(); 

  // Get candidate Ptms with similar mass shifts
  PtmPtrVec match_ptm_ptr_vec = local_util::getPtmPtrVecByMass(unexp_shift_mass, err_tole, ptm_ptr_vec_);

  // if there is a match, find the best ptm and its best similarity score
  if (match_ptm_ptr_vec.size() != 0) {
    // compute similarity score for each possible site of the PTM
    LocalResultPtr local_result_ptr = compOnePtmScr(form_ptr, extend_ms_ptr_vec, unexp_shift_mass, match_ptm_ptr_vec);
    return local_result_ptr;
  }
  return nullptr;
}

ProteoformPtr LocalProcessor::createProteoformPtr(ProteoformPtr base_form_ptr, 
                                                  double shift_mass, LocalResultPtr result_ptr) {
  if (result_ptr->match_score_ <= 0) {
    return nullptr;
  }
  int bgn, end;
  double conf;
  local_util::scrFilter(result_ptr->scr_vec_, bgn, end, conf, mng_ptr_->threshold_);

  if (bgn == -1) return nullptr;

  LocalAnnoPtr anno = std::make_shared<LocalAnno>(bgn, end, conf, result_ptr->scr_vec_,
                                                  result_ptr->match_score_, result_ptr->ptm_ptr_);

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

LocalResultPtr LocalProcessor::compOnePtmScr(ProteoformPtr base_form_ptr, 
                                             const ExtendMsPtrVec & extend_ms_ptr_vec,
                                             double shift_mass, 
                                             PtmPtrVec & ptm_ptr_vec) {
  // generate a new shift using the sum of shift masses
  int left_pos = 0;
  int right_pos = 1;
  MassShiftPtr unexp_shift_ptr = local_util::geneMassShift(left_pos, right_pos, shift_mass,
                                                           AlterType::UNEXPECTED);
  MassShiftPtrVec all_shift_ptr_vec = base_form_ptr->getMassShiftPtrVec();
  all_shift_ptr_vec.push_back(unexp_shift_ptr);

  // create a proteoform with all modifications, including the unexpected mass
  // shift
  ProteoformPtr shift_form_ptr
    = std::make_shared<Proteoform>(base_form_ptr->getFastaSeqPtr(),
                                   base_form_ptr->getProtModPtr(),
                                   base_form_ptr->getStartPos(),
                                   base_form_ptr->getEndPos(),
                                   base_form_ptr->getResSeqPtr(),
                                   all_shift_ptr_vec);

  int len = shift_form_ptr->getLen();

  LocalResultPtr result_ptr = std::make_shared<LocalResult>();
  // current score vector
  std::vector<double> cur_scr_vec;

  // check candidate PTMs one by one
  for (size_t i = 0; i < ptm_ptr_vec.size(); i++) {
    cur_scr_vec.clear();
    // number of modification sites
    int cur_count = 0;
    // best best match
    int cur_best_match = -1;
    for (int j = 0; j < len; j++) {
      if (modifiable(shift_form_ptr, j, ptm_ptr_vec[i])) {
        cur_count++;
        unexp_shift_ptr->setLeftBpPos(j);
        unexp_shift_ptr->setRightBpPos(j + 1);
        int match = static_cast<int>(local_util::compMatchFragNum(shift_form_ptr,
                                                                  extend_ms_ptr_vec,
                                                                  mng_ptr_->min_mass_));

        cur_scr_vec.push_back(std::pow(mng_ptr_->p1_, len - match) * std::pow(mng_ptr_->p2_, match));
        if (match > cur_best_match) {
          cur_best_match = match;
        }
      } else {
        cur_scr_vec.push_back(0.0);
      }
    }

    // If the best match of the PTM has a higher score than that of previous
    // PTMs, it is reported as the best PTM. 
    if (cur_count > 0 && cur_best_match > result_ptr->match_score_) {
      result_ptr->match_score_ = cur_best_match;
      result_ptr->ptm_ptr_ = ptm_ptr_vec[i];
      result_ptr->scr_vec_ = cur_scr_vec;
    }
  }
  
  // normalize
  if (result_ptr->match_score_ > 0) {
    local_util::normalize(result_ptr->scr_vec_);
    ProteoformPtr local_form_ptr = createProteoformPtr(base_form_ptr, shift_mass, result_ptr);
    if (local_form_ptr != nullptr) {
      result_ptr->form_ptr_ = local_form_ptr;
      return result_ptr;
    }
  }

  return nullptr;
}

bool LocalProcessor::modifiable(ProteoformPtr proteoform_ptr, int i, PtmPtr ptm_ptr) {
  if (ptm_ptr == nullptr) return true;

  MassShiftPtrVec fixed_shift_vec = proteoform_ptr->getMassShiftPtrVec(AlterType::FIXED);

  for (size_t k = 0; k < fixed_shift_vec.size(); k++) {
    if (fixed_shift_vec[k]->getLeftBpPos() <= i && i < fixed_shift_vec[k]->getRightBpPos()) {
      return false;
    }
  }

  int start = proteoform_ptr->getStartPos();
  int end = proteoform_ptr->getEndPos();

  ResiduePtr residue_ptr = proteoform_ptr->getResSeqPtr()->getResiduePtr(i);

  ModPtrVec mod_list;

  if (i == 0 && start == 0) {
    mod_list = mod_list_N_;
  } else if (i + start == end) {
    mod_list = mod_list_C_;
  } else {
    mod_list = mod_list_any_;
  }

  for (size_t j = 0; j < mod_list.size(); j++) {
    if (mod_list[j]->getOriResiduePtr()->isSame(residue_ptr) &&
        mod_list[j]->getModResiduePtr()->getPtmPtr()->isSame(ptm_ptr))
      return true;
  }

  return false;
}

/*
void findBestPtm(FastaSeqPtr fasta_seq_ptr, ProtModPtr prot_mod_ptr, 
                 int start_pos, int end_pos, 
                 ResSeqPtr res_seq_ptr, 
                 MassShiftPtrVec &expected_shift_vec, 
                 double unexpected_mass_shift, 
                 PtmPtrVec &ptm_ptr_vec) {

  if (ptm_ptr_vec.size() == 0) {
    LOG_ERROR("Ptm list is empty!");
    exit(EXIT_FAILURE);
  }

  MassShiftPtrVec all_shift_vec = expected_shift_vec;
  all_shift_vec.push_back(unexpe_mass_shift_ptr);
  std::sort(all_shift_vec.begin(), all_shift_vec.end(), MassShift::cmpPosInc);
  ProteoformPtr one_shift_proteoform
      = std::make_shared<Proteoform>(fast_seq_ptr, 
                                     prot_mod_ptr, 
                                     start_pos, 
                                     end_pos, 
                                     res_seq_ptr, 
                                     all_shift_vec);

  // confidence score
  std::vector<double> scr_vec;
  int len = end_pos - start_pos + 1;
  for (size_t i = 0; i < ptm_ptr_vec.size(); i++) {
    scr_vec.clear();
    for (int j = 0; j < len; j++) {
      if (modifiable(proteoform, j, ptm_vec[i])) {
        mass_shift->setLeftBpPos(j);
        mass_shift->setRightBpPos(j + 1);
        int match = static_cast<int>(local_util::compMatchFragNum(proteoform,
                                                                  extend_ms_ptr_vec,
                                                                  mng_ptr_->min_mass_));

        scr_vec.push_back(std::pow(mng_ptr_->p1_, n - match) * std::pow(mng_ptr_->p2_, match));
      } else {
        scr_vec.push_back(0.0);
      }
    }
    scr_vec2d.push_back(scr_vec);
    temp.push_back(std::accumulate(scr_vec.begin(), scr_vec.end(), 0.0));
  }


}
*/

/*
ProteoformPtr LocalProcessor::twoPtmTermAdjust(ProteoformPtr proteoform, int num_match,
                                               const ExtendMsPtrVec & extend_ms_ptr_vec, double prec_mass,
                                               double & mass1, double & mass2) {
  int n_trunc_min, n_trunc_max, c_trunc_min, c_trunc_max;
  getNtermTruncRange(proteoform, extend_ms_ptr_vec, n_trunc_min, n_trunc_max);
  getCtermTruncRange(proteoform, extend_ms_ptr_vec, c_trunc_min, c_trunc_max);

  MassShiftPtrVec ori_mass_shift_vec = proteoform->getMassShiftPtrVec();

  MassShiftPtrVec ori_fix_mass_shift_vec
      = local_util::massShiftFilter(ori_mass_shift_vec, AlterType::UNEXPECTED);

  MassShiftPtr mass_shift1 = proteoform->getMassShiftPtrVec(AlterType::UNEXPECTED)[0];
  MassShiftPtr mass_shift2 = proteoform->getMassShiftPtrVec(AlterType::UNEXPECTED)[1];

  double err = prec_mass * mng_ptr_->ppo_;

  int ori_start = proteoform->getStartPos();
  int ori_end = proteoform->getEndPos();

  double ori_mass1 = mass_shift1->getMassShift();
  double ori_mass2 = mass_shift2->getMassShift();
  mass1 = ori_mass1;
  mass2 = ori_mass2;

  std::vector<bool> ptm1_known_vec, ptm2_known_vec;
  std::vector<int> c_vec, n_vec;
  std::string n_seq, c_seq;
  std::vector<double> raw_scr_vec;

  for (int i = n_trunc_min; i <= n_trunc_max; i++) {
    for (int j = c_trunc_min; j <= c_trunc_max; j++) {
      if (i < 0) {
        n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start + i, -i);
        mass1 = ori_mass1 - residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      } else {
        n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start, i);
        mass1 = ori_mass1 + residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }
      if (j >= 0) {
        c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1, j);
        mass2 = ori_mass2 - residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      } else {
        c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1 + j, -j);
        mass2 = ori_mass2 + residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }

      if ((mass1 > mng_ptr_->max_ptm_mass_ || mass1 < mng_ptr_->min_ptm_mass_
           || mass2 > mng_ptr_->max_ptm_mass_ || mass2 < mng_ptr_->min_ptm_mass_) && (i != 0 || j != 0))
        continue;

      n_vec.push_back(i);
      c_vec.push_back(j);

      PtmPairVec ptm_pair_vec = local_util::getPtmPairVecByMass(mass1, mass2, err, ptm_pair_vec_);
      if (ptm_pair_vec.size() == 0) {
        ptm_pair_vec.push_back(std::make_pair(nullptr, nullptr));
      }

      ptm1_known_vec.push_back(local_util::getPtmPtrVecByMass(mass1, err, ptm_vec_).size() > 0);
      ptm2_known_vec.push_back(local_util::getPtmPtrVecByMass(mass2, err, ptm_vec_).size() > 0);

      mass_shift1 = local_util::geneMassShift(mass_shift1, mass1, AlterType::UNEXPECTED);
      mass_shift2 = local_util::geneMassShift(mass_shift2, mass2, AlterType::UNEXPECTED);

      MassShiftPtrVec new_mass_shift_vec = local_util::copyMassShiftVec(ori_fix_mass_shift_vec);
      new_mass_shift_vec.push_back(mass_shift1);
      new_mass_shift_vec.push_back(mass_shift2);

      std::sort(new_mass_shift_vec.begin(), new_mass_shift_vec.end(), MassShift::cmpPosInc);

      ProteoformPtr tmp_proteoform = proteoform_factory::geneProteoform(proteoform,
                                                                        ori_start + i,
                                                                        ori_end + j,
                                                                        new_mass_shift_vec,
                                                                        mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      double raw_scr;
      compTwoPtmScr(tmp_proteoform, num_match, extend_ms_ptr_vec, prec_mass, raw_scr, ptm_pair_vec);
      raw_scr_vec.push_back(raw_scr);
    }
  }

  for (size_t i = 0; i < raw_scr_vec.size(); i++) {
    if (ptm1_known_vec[i] && ptm2_known_vec[i]) {
      raw_scr_vec[i] = raw_scr_vec[i] * mng_ptr_->theta_ * mng_ptr_->theta_;
    } else if (ptm1_known_vec[i] || ptm2_known_vec[i]) {
      raw_scr_vec[i] = raw_scr_vec[i] * mng_ptr_->theta_ * (1 - mng_ptr_->theta_);
    } else {
      raw_scr_vec[i] = raw_scr_vec[i] * (1 - mng_ptr_->theta_) * (1 - mng_ptr_->theta_);
    }
  }

  int idx = std::distance(raw_scr_vec.begin(), std::max_element(raw_scr_vec.begin(), raw_scr_vec.end()));

  if (n_vec[idx] < 0) {
    n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start + n_vec[idx], -n_vec[idx]);
    mass1 = ori_mass1 - residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  } else {
    n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start, n_vec[idx]);
    mass1 = ori_mass1 + residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  }

  if (c_vec[idx] >= 0) {
    c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1, c_vec[idx]);
    mass2 = ori_mass2 - residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  } else {
    c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1 + c_vec[idx], -c_vec[idx]);
    mass2 = ori_mass2 + residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  }

  mass_shift1 = local_util::geneMassShift(mass_shift1, mass1, AlterType::UNEXPECTED);
  mass_shift2 = local_util::geneMassShift(mass_shift2, mass2, AlterType::UNEXPECTED);

  MassShiftPtrVec new_mass_shift_vec = local_util::copyMassShiftVec(ori_fix_mass_shift_vec);
  new_mass_shift_vec.push_back(mass_shift1);
  new_mass_shift_vec.push_back(mass_shift2);

  std::sort(new_mass_shift_vec.begin(), new_mass_shift_vec.end(), MassShift::cmpPosInc);

  return proteoform_factory::geneProteoform(proteoform,
                                            ori_start + n_vec[idx],
                                            ori_end + c_vec[idx],
                                            new_mass_shift_vec,
                                            mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
}


// similar to processOneKnownPtm, we might get a nullptr from this function
ProteoformPtr LocalProcessor::processTwoKnownPtm(PrsmPtr prsm) {
  ProteoformPtr ori_prot_form = prsm->getProteoformPtr();

  MassShiftPtrVec ori_mass_shift_vec = ori_prot_form->getMassShiftPtrVec();

  MassShiftPtrVec unexpected_shift_vec
      = ori_prot_form->getMassShiftPtrVec(AlterType::UNEXPECTED);

  MassShiftPtrVec expected_shift_vec = local_util::massShiftFilter(ori_mass_shift_vec,
                                                                   AlterType::UNEXPECTED);

  MassShiftPtr unexpected_shift1 = local_util::geneMassShift(unexpected_shift_vec[0],
                                                             unexpected_shift_vec[0]->getMassShift(),
                                                             AlterType::UNEXPECTED);

  MassShiftPtr unexpected_shift2;

  if (unexpected_shift_vec.size() == 2) {
    unexpected_shift2 = local_util::geneMassShift(unexpected_shift_vec[1],
                                                  unexpected_shift_vec[1]->getMassShift(),
                                                  AlterType::UNEXPECTED);
  } else {
    unexpected_shift2 = local_util::geneMassShift(unexpected_shift_vec[0],
                                                  0.0, AlterType::UNEXPECTED);
  }

  double shift_mass1 = unexpected_shift1->getMassShift();

  double shift_mass2 = unexpected_shift2->getMassShift();

  expected_shift_vec.push_back(unexpected_shift1);

  expected_shift_vec.push_back(unexpected_shift2);

  std::sort(expected_shift_vec.begin(), expected_shift_vec.end(), MassShift::cmpPosInc);

  ProteoformPtr two_shift_proteoform
      = std::make_shared<Proteoform>(ori_prot_form->getFastaSeqPtr(),
                                     ori_prot_form->getProtModPtr(),
                                     ori_prot_form->getStartPos(),
                                     ori_prot_form->getEndPos(),
                                     ori_prot_form->getResSeqPtr(),
                                     expected_shift_vec);

  double err = prsm->getAdjustedPrecMass() * mng_ptr_->ppo_;

  PtmPairVec ptm_pair_vec = local_util::getPtmPairVecByMass(shift_mass1, shift_mass2,
                                                            err, ptm_pair_vec_);

  if (ptm_pair_vec.size() == 0) {
    if (two_shift_proteoform->getProtModPtr()->getType() == ProtModBase::getType_NME_ACETYLATION()
        || two_shift_proteoform->getProtModPtr()->getType() == ProtModBase::getType_M_ACETYLATION()) {
      MassShiftPtrVec expected_shift_vec = local_util::massShiftFilter(ori_mass_shift_vec,
                                                                       AlterType::UNEXPECTED);
      for (size_t k = 0; k < expected_shift_vec.size(); k++) {
        if (expected_shift_vec[k]->getTypePtr() == AlterType::PROTEIN_VARIABLE) {
          shift_mass1 += expected_shift_vec[k]->getMassShift();
          expected_shift_vec.erase(expected_shift_vec.begin() + k);
          break;
        }
      }

      unexpected_shift1 = local_util::geneMassShift(unexpected_shift1, shift_mass1,
                                                    AlterType::UNEXPECTED);

      expected_shift_vec.push_back(unexpected_shift1);

      expected_shift_vec.push_back(unexpected_shift2);

      std::sort(expected_shift_vec.begin(), expected_shift_vec.end(), MassShift::cmpPosInc);

      two_shift_proteoform = std::make_shared<Proteoform>(ori_prot_form->getFastaSeqPtr(),
                                                          ori_prot_form->getProtModPtr(),
                                                          ori_prot_form->getStartPos(),
                                                          ori_prot_form->getEndPos(),
                                                          ori_prot_form->getResSeqPtr(),
                                                          expected_shift_vec);
    }

    twoPtmTermAdjust(two_shift_proteoform, prsm->getMatchPeakNum(), prsm->getRefineMsPtrVec(),
                     prsm->getAdjustedPrecMass(), shift_mass1, shift_mass2);

    ptm_pair_vec = local_util::getPtmPairVecByMass(shift_mass1, shift_mass2, err, ptm_pair_vec_);
  }

  if (ptm_pair_vec.size() == 0) return nullptr;

  double raw_scr;
  ExtendMsPtrVec extend_ms_ptr_vec = prsm->getRefineMsPtrVec();

  compTwoPtmScr(two_shift_proteoform, prsm->getMatchPeakNum(), extend_ms_ptr_vec,
                prsm->getAdjustedPrecMass(),
                raw_scr, ptm_pair_vec);

  // can be explained by a variable ptm, but no modifiable site
  if (raw_scr == 0) return nullptr;

  two_shift_proteoform->setProteoClusterId(ori_prot_form->getProteoClusterId());
  two_shift_proteoform->setProtId(ori_prot_form->getProtId());

  PtmPtr ptm1 = ptm_pair_vec[0].first;
  PtmPtr ptm2 = ptm_pair_vec[0].second;
  raw_scr = raw_scr * mng_ptr_->theta_ * mng_ptr_->theta_ * (1 - mng_ptr_->beta_);

  std::vector<double> empty_scr_vec;
  LocalAnnoPtr anno1 = std::make_shared<LocalAnno>(0, 0, 0, empty_scr_vec, raw_scr, ptm1);
  two_shift_proteoform->getMassShiftPtrVec(AlterType::UNEXPECTED)[0]->getAlterPtr(0)->setLocalAnno(anno1);

  LocalAnnoPtr anno2 = std::make_shared<LocalAnno>(0, 0, 0, empty_scr_vec, raw_scr, ptm2);
  two_shift_proteoform->getMassShiftPtrVec(AlterType::UNEXPECTED)[1]->getAlterPtr(0)->setLocalAnno(anno2);

  return compSplitPoint(two_shift_proteoform, prsm->getMatchPeakNum(),
                        extend_ms_ptr_vec, prsm->getAdjustedPrecMass());
}

PrsmPtr LocalProcessor::processTwoPtm(PrsmPtr prsm) {
  int ori_num_match_ion = local_util::compMatchFragNum(prsm->getProteoformPtr(),
                                                       prsm->getRefineMsPtrVec(),
                                                       mng_ptr_->min_mass_);

  ProteoformPtr two_known_prsm = processTwoKnownPtm(prsm);

  if (two_known_prsm != nullptr) {
    int new_num_match_ion = local_util::compMatchFragNum(two_known_prsm,
                                                         prsm->getRefineMsPtrVec(),
                                                         mng_ptr_->min_mass_);
    if (new_num_match_ion > ori_num_match_ion - mng_ptr_->DESC_MATCH_LIMIT_
        && new_num_match_ion > ori_num_match_ion * mng_ptr_->desc_ratio_) {
      prsm->setProteoformPtr(two_known_prsm, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      return prsm;
    }
  }

  ProteoformPtr one_known_prsm = processOneKnownPtm(prsm);

  if (one_known_prsm != nullptr) {
    double new_num_match_ion = local_util::compMatchFragNum(one_known_prsm,
                                                            prsm->getRefineMsPtrVec(),
                                                            mng_ptr_->min_mass_);
    if (new_num_match_ion > ori_num_match_ion - mng_ptr_->DESC_MATCH_LIMIT_
        && new_num_match_ion > ori_num_match_ion * mng_ptr_->desc_ratio_) {
      prsm->setProteoformPtr(one_known_prsm, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      return prsm;
    }
  }

  return prsm;
}


void LocalProcessor::compTwoPtmScr(ProteoformPtr proteoform, int num_match,
                                   const ExtendMsPtrVec & extend_ms_ptr_vec, double prec_mass,
                                   double & raw_scr, PtmPairVec & ptm_pair_vec) {
  MassShiftPtrVec ori_mass_shift_vec = proteoform->getMassShiftPtrVec();
  double shift_mass1 = proteoform->getMassShiftPtrVec(AlterType::UNEXPECTED)[0]->getMassShift();
  double shift_mass2 = proteoform->getMassShiftPtrVec(AlterType::UNEXPECTED)[1]->getMassShift();
  ProteoformPtr no_shift_proteoform
      = proteoform_factory::geneProteoform(proteoform,
                                           proteoform->getStartPos(),
                                           proteoform->getEndPos(),
                                           local_util::massShiftFilter(ori_mass_shift_vec, AlterType::UNEXPECTED),
                                           mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  std::vector<double> scr_vec;

  for (size_t k = 0; k < ptm_pair_vec.size(); k++) {
    local_util::ptmMassAdjust(shift_mass1, shift_mass2, ptm_pair_vec[k].first, ptm_pair_vec[k].second);
    scr_vec.push_back(dpTwoPtmScr(no_shift_proteoform, num_match, extend_ms_ptr_vec, prec_mass,
                                  shift_mass1, shift_mass2, ptm_pair_vec[k].first, ptm_pair_vec[k].second));
  }

  int idx = std::distance(scr_vec.begin(), std::max_element(scr_vec.begin(), scr_vec.end()));
  raw_scr = scr_vec[idx];
  PtmPtr p1 = ptm_pair_vec[idx].first;
  PtmPtr p2 = ptm_pair_vec[idx].second;
  ptm_pair_vec.clear();
  ptm_pair_vec.push_back(std::make_pair(p1, p2));
}

double LocalProcessor::dpTwoPtmScr(ProteoformPtr proteoform, int h,
                                   const ExtendMsPtrVec & extend_ms_ptr_vec,
                                   double prec_mass, double mass1, double mass2, PtmPtr ptm1, PtmPtr ptm2) {
  int g = proteoform->getLen();
  if (g <= 0) {
    return 0;
  }
  double scr = 0.0;
  int count = 0;

  for (size_t k = 0; k < extend_ms_ptr_vec.size(); k++) {
    std::vector<std::vector<double>> b_table(3);
    b_table[0] = local_util::geneNTheoMass(proteoform, extend_ms_ptr_vec[k], mng_ptr_->min_mass_);
    local_util::fillTableB(b_table, mass1, mass2);

    std::vector<std::vector<int>> s_table(3);
    local_util::fillTableS(b_table, s_table, extend_ms_ptr_vec[k], prec_mass);

    // fill D(f,g,h)
    int d_table[3][g + 1][h + 1];
    memset(d_table, 0, sizeof(int) * 3 * (g + 1) * (h + 1));
    d_table[0][0][0] = 1;

    for (int i = 1; i <= g; i++) {
      for (int j = 0; j <= h; j++) {
        if (j >= s_table[0][i - 1]) {
          d_table[0][i][j] = d_table[0][i - 1][j - s_table[0][i - 1]];
        } else {
          d_table[0][i][j] = 0;
        }
      }
    }

    for (int i = 1; i <= g; i++) {
      for (int j = 0; j <= h; j++) {
        if (modifiable(proteoform, i - 1, ptm1) && j >= s_table[1][i - 1]) {
          d_table[1][i][j] = d_table[0][i - 1][j - s_table[1][i - 1]] + d_table[1][i - 1][j - s_table[1][i - 1]];
        } else if (j >= s_table[1][i - 1]) {
          d_table[1][i][j] = d_table[1][i - 1][j - s_table[1][i - 1]];
        } else {
          d_table[1][i][j] = 0;
        }
      }
    }

    for (int i = 1; i <= g; i++) {
      for (int j = 0; j <= h; j++) {
        if (modifiable(proteoform, i - 1, ptm2) && j >= s_table[2][i - 1]) {
          d_table[2][i][j] = d_table[1][i - 1][j - s_table[2][i - 1]] + d_table[2][i - 1][j - s_table[2][i - 1]];
        } else if (j >= s_table[2][i - 1]) {
          d_table[2][i][j] = d_table[2][i - 1][j - s_table[2][i - 1]];
        } else {
          d_table[2][i][j] = 0;
        }
      }
    }

    for (int i = 0; i <= h; i++) {
      count += d_table[2][g][i];
      scr += d_table[2][g][i] * std::pow(mng_ptr_->p1_, g - i) * std::pow(mng_ptr_->p2_, i);
    }
  }

  return scr / count;
}

ProteoformPtr LocalProcessor::compSplitPoint(ProteoformPtr proteoform, int h,
                                             const ExtendMsPtrVec & extend_ms_ptr_vec,
                                             double prec_mass) {
  MassShiftPtrVec ori_mass_shift_vec = proteoform->getMassShiftPtrVec();

  MassShiftPtr shift_ptr1 = proteoform->getMassShiftPtrVec(AlterType::UNEXPECTED)[0];
  MassShiftPtr shift_ptr2 = proteoform->getMassShiftPtrVec(AlterType::UNEXPECTED)[1];

  double mass1 = shift_ptr1->getMassShift();
  double mass2 = shift_ptr2->getMassShift();

  PtmPtr ptm1 = shift_ptr1->getAlterPtr(0)->getLocalAnno()->getPtmPtr();
  PtmPtr ptm2 = shift_ptr2->getAlterPtr(0)->getLocalAnno()->getPtmPtr();

  local_util::ptmMassAdjust(mass1, mass2, ptm1, ptm2);

  shift_ptr1 = local_util::geneMassShift(shift_ptr1, mass1, AlterType::UNEXPECTED);

  shift_ptr2 = local_util::geneMassShift(shift_ptr2, mass2, AlterType::UNEXPECTED);

  MassShiftPtrVec expected_shift_vec = local_util::massShiftFilter(ori_mass_shift_vec, AlterType::UNEXPECTED);

  expected_shift_vec.push_back(shift_ptr1);

  expected_shift_vec.push_back(shift_ptr2);

  int prot_cluster_id = proteoform->getProteoClusterId();

  int prot_id = proteoform->getProtId();

  proteoform = proteoform_factory::geneProteoform(proteoform,
                                                  proteoform->getStartPos(),
                                                  proteoform->getEndPos(),
                                                  expected_shift_vec,
                                                  mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  proteoform->setProteoClusterId(prot_cluster_id);

  proteoform->setProtId(prot_id);

  int g = proteoform->getLen();

  ProteoformPtr no_shift_proteoform
      = proteoform_factory::geneProteoform(proteoform,
                                           proteoform->getStartPos(),
                                           proteoform->getEndPos(),
                                           local_util::massShiftFilter(ori_mass_shift_vec, AlterType::UNEXPECTED),
                                           mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  std::vector<double> split_scr_vec;
  for (int k = 1; k < g; k++) {
    double scr = 0.0;
    for (size_t i = 0; i < extend_ms_ptr_vec.size(); i++) {
      std::vector<std::vector<double>> b_table(3);
      b_table[0] = local_util::geneNTheoMass(no_shift_proteoform, extend_ms_ptr_vec[i], mng_ptr_->min_mass_);
      local_util::fillTableB(b_table, mass1, mass2);

      std::vector<std::vector<int>> s_table(3);
      local_util::fillTableS(b_table, s_table, extend_ms_ptr_vec[i], prec_mass);

      int d_table[3][g + 1][h + 1];

      memset(d_table, 0, sizeof(int) * 3 * (g + 1) * (h + 1));
      d_table[0][0][0] = 1;

      for (int i = 1; i <= g; i++) {
        for (int j = 0; j <= h; j++) {
          if (j >= s_table[0][i - 1]) {
            d_table[0][i][j] = d_table[0][i - 1][j - s_table[0][i - 1]];
          } else {
            d_table[0][i][j] = 0;
          }
        }
      }

      for (int i = 1; i <= g; i++) {
        for (int j = 0; j <= h; j++) {
          if (modifiable(proteoform, i - 1, ptm1) && j >= s_table[1][i - 1] && i <= k) {
            d_table[1][i][j] = d_table[0][i - 1][j - s_table[1][i - 1]] + d_table[1][i - 1][j - s_table[1][i - 1]];
          } else if (j >= s_table[1][i - 1]) {
            d_table[1][i][j] = d_table[1][i - 1][j - s_table[1][i - 1]];
          } else {
            d_table[1][i][j] = 0;
          }
        }
      }

      for (int i = 1; i <= g; i++) {
        for (int j = 0; j <= h; j++) {
          if (modifiable(proteoform, i - 1, ptm2) && j >= s_table[2][i - 1] && i > k) {
            d_table[2][i][j] = d_table[1][i - 1][j - s_table[2][i - 1]] + d_table[2][i - 1][j - s_table[2][i - 1]];
          } else if (j >= s_table[2][i - 1]) {
            d_table[2][i][j] = d_table[2][i - 1][j - s_table[2][i - 1]];
          } else {
            d_table[2][i][j] = 0;
          }
        }
      }

      for (int i = 0; i <= h; i++) {
        scr += d_table[2][g][i] * std::pow(mng_ptr_->p1_, i) * std::pow(mng_ptr_->p2_, i);
      }
    }
    split_scr_vec.push_back(scr);
  }

  int split_point = std::distance(split_scr_vec.begin(), std::max_element(split_scr_vec.begin(), split_scr_vec.end()));

  double split_max = *std::max_element(split_scr_vec.begin(), split_scr_vec.end());

  double split_scr = 0.0;

  int split_end = split_scr_vec.size();

  for (; split_end > 0; split_end--) {
    if (split_scr_vec[split_end] == split_max)
      break;
  }

  for (int i = split_point; i <= split_end; i++) {
    split_scr += split_scr_vec[i];
  }

  split_point = (split_point + split_end) / 2;

  split_scr = split_scr / std::accumulate(split_scr_vec.begin(), split_scr_vec.end(), 0.0);

  if (split_scr <= mng_ptr_->threshold_) {
    return nullptr;
  }

  std::vector<double> ptm_scr;

  for (int i = 0; i <= split_point; i++) {
    shift_ptr1->setLeftBpPos(i);
    shift_ptr1->setRightBpPos(i + 1);
    if (modifiable(proteoform, i, ptm1)) {
      int match = static_cast<int>(local_util::compMatchFragNum(proteoform,
                                                                extend_ms_ptr_vec,
                                                                mng_ptr_->min_mass_));
      ptm_scr.push_back(std::pow(mng_ptr_->p1_, split_point - match) * std::pow(mng_ptr_->p2_, match));
    } else {
      ptm_scr.push_back(0.0);
    }
  }

  local_util::normalize(ptm_scr);
  int bgn, end;
  double conf;
  std::transform(ptm_scr.begin(), ptm_scr.end(), ptm_scr.begin(),
                 std::bind1st(std::multiplies<double>(), split_scr));
  local_util::scrFilter(ptm_scr, bgn, end, conf, mng_ptr_->threshold_);

  if (bgn == -1) {
    return nullptr;
  } else {
    LocalAnnoPtr anno1 = std::make_shared<LocalAnno>(bgn, end, conf, ptm_scr, 0, ptm1);
    shift_ptr1->getAlterPtr(0)->setLocalAnno(anno1);
    shift_ptr1->setLeftBpPos(anno1->getLeftBpPos());
    shift_ptr1->setRightBpPos(anno1->getRightBpPos() + 1);
  }

  ptm_scr.clear();
  int len = proteoform->getLen();
  for (int i = split_point + 1; i < len; i++) {
    shift_ptr2->setLeftBpPos(i);
    shift_ptr2->setRightBpPos(i + 1);
    if (modifiable(proteoform, i, ptm2)) {
      int match = static_cast<int>(local_util::compMatchFragNum(proteoform,
                                                                extend_ms_ptr_vec,
                                                                mng_ptr_->min_mass_));
      ptm_scr.push_back(std::pow(mng_ptr_->p1_, len - match) * std::pow(mng_ptr_->p2_, match));
    } else {
      ptm_scr.push_back(0.0);
    }
  }

  local_util::normalize(ptm_scr);

  std::transform(ptm_scr.begin(), ptm_scr.end(), ptm_scr.begin(),
                 std::bind1st(std::multiplies<double>(), split_scr));

  local_util::scrFilter(ptm_scr, bgn, end, conf, mng_ptr_->threshold_);

  if (bgn == -1) {
    return nullptr;
  } else {
    LocalAnnoPtr anno2 = std::make_shared<LocalAnno>(split_point + bgn + 1, split_point + end + 1, conf, ptm_scr, 0, ptm2);
    shift_ptr2->getAlterPtr(0)->setLocalAnno(anno2);
    shift_ptr2->setLeftBpPos(anno2->getLeftBpPos());
    shift_ptr2->setRightBpPos(anno2->getRightBpPos() + 1);
  }

  return proteoform;
}
*/

}  // namespace toppic

