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
#include <iomanip>
#include <limits>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/base/ptm_util.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/residue_util.hpp"

#include "ms/spec/msalign_util.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"

#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"

#include "stat/local/local_util.hpp"
#include "stat/local/local_proteoform.hpp"
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
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileNameWithFolder();

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

        if (prsm_ptr->getProteoformPtr()->getAlterNum(AlterType::UNEXPECTED) > 0) {
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
  int mass_shift_num = prsm->getProteoformPtr()->getAlterNum(AlterType::UNEXPECTED);
  double err_tole = mng_ptr_->peak_tole_ptr_->compStrictErrorTole(prsm->getAdjustedPrecMass());
  if (mass_shift_num == 1) {
    double mass = prsm->getProteoformPtr()->getMassShiftPtrVec(AlterType::UNEXPECTED)[0]->getMassShift();
    // if the mass shift is between [-1, 1], we don't characterize it
    if (std::abs(mass) <= 1 + err_tole) return prsm;
    //LOG_DEBUG("start one mass shift analysis");
    return processOneMassShift(prsm);
  } else if (mass_shift_num == 2) {
    //LOG_DEBUG("start two mass shifts analysis");
    return processTwoMassShifts(prsm);
  }
  return prsm;
}

PrsmPtr LocalProcessor::processOneMassShift(PrsmPtr prsm) {
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
      prsm->setAdjustedPrecMass(one_known_proteoform->getMass());
      return prsm;
    }
  }

  /*
  TWO PTM explanation may introduce some problems of proteoform annotation.
  We will add the function after proteoform annotation is reviewed. 
  // Check if the mass shift can be explained by two common PTMs
  LOG_DEBUG("start two known ptm localization")
  ProteoformPtr two_known_proteoform = processTwoKnownPtms(prsm);

  if (two_known_proteoform != nullptr) {
    double new_num_match_ion = local_util::compMatchFragNum(two_known_proteoform,
                                                            prsm->getRefineMsPtrVec(),
                                                            mng_ptr_->min_mass_);
    //LOG_DEBUG("orig score " << ori_num_match_ion << " new match num " << new_num_match_ion);
    //LOG_DEBUG("two known form " << two_known_proteoform->getProteinMatchSeq());
    if (new_num_match_ion > ori_num_match_ion - mng_ptr_->DESC_MATCH_LIMIT_
        && new_num_match_ion > ori_num_match_ion * mng_ptr_->desc_ratio_) {
      two_known_proteoform->setProteoClusterId(prsm->getProteoformPtr()->getProteoClusterId());
      two_known_proteoform->setProtId(prsm->getProteoformPtr()->getProtId());
      prsm->setProteoformPtr(two_known_proteoform, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      return prsm;
    }
  }
  */
  return prsm;
}

PrsmPtr LocalProcessor::processTwoMassShifts(PrsmPtr prsm) {
  int ori_num_match_ion = local_util::compMatchFragNum(prsm->getProteoformPtr(),
                                                       prsm->getRefineMsPtrVec(),
                                                       mng_ptr_->min_mass_);

  ProteoformPtr two_known_proteoform = processTwoKnownPtms(prsm);

  if (two_known_proteoform != nullptr) {
    int new_num_match_ion = local_util::compMatchFragNum(two_known_proteoform,
                                                         prsm->getRefineMsPtrVec(),
                                                         mng_ptr_->min_mass_);
    if (new_num_match_ion > ori_num_match_ion - mng_ptr_->DESC_MATCH_LIMIT_
        && new_num_match_ion > ori_num_match_ion * mng_ptr_->desc_ratio_) {
      two_known_proteoform->setProteoClusterId(prsm->getProteoformPtr()->getProteoClusterId());
      two_known_proteoform->setProtId(prsm->getProteoformPtr()->getProtId());
      prsm->setProteoformPtr(two_known_proteoform, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      prsm->setAdjustedPrecMass(two_known_proteoform->getMass());
      return prsm;
    }
  }

  return prsm;
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

  if (i == 0) {
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


// we will get a nullptr if the mass shift can't be explained by a variable ptm
ProteoformPtr LocalProcessor::processOneKnownPtm(PrsmPtr prsm_ptr) {
  //get canidate forms
  ProteoformPtr ori_form_ptr = prsm_ptr->getProteoformPtr();
  ProteoformPtrVec cand_form_vec = local_proteoform::getAllCandidateForms(ori_form_ptr, mng_ptr_);

  ExtendMsPtrVec extend_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  double adjust_prec_mass = prsm_ptr->getAdjustedPrecMass();
  double err_tole = mng_ptr_->peak_tole_ptr_->compStrictErrorTole(adjust_prec_mass);

  int best_match_score = 0;
  ProteoformPtr best_form_ptr = nullptr;
  PtmPtr best_ptm_ptr = nullptr;

  //LOG_DEBUG("form number " << cand_form_vec.size());
  for (size_t i = 0; i < cand_form_vec.size(); i++) {
    ProteoformPtr cand_form_ptr = cand_form_vec[i];
    double unexp_shift_mass = adjust_prec_mass - cand_form_ptr->getMass(); 
    LOG_DEBUG(std::setprecision(10) << "adjust prec mass " << adjust_prec_mass << " form mass " << cand_form_ptr->getMass());
    // Get candidate Ptms with similar mass shifts
    PtmPtrVec match_ptm_ptr_vec = local_util::getPtmPtrVecByMass(unexp_shift_mass, err_tole, ptm_ptr_vec_);
    LOG_DEBUG("ptm number " << match_ptm_ptr_vec.size() << " shift " << unexp_shift_mass);

    // if there is a match, find the best ptm and its best similarity score
    for (size_t j = 0; j < match_ptm_ptr_vec.size(); j++) {
      PtmPtr ptm_ptr = match_ptm_ptr_vec[j];
      // compute similarity score for each possible site of the PTM
      int score = compOnePtmScr(cand_form_ptr, extend_ms_ptr_vec, ptm_ptr); 
      //LOG_DEBUG("score " << score);
      if (score > best_match_score) {
        best_match_score = score;
        best_form_ptr = cand_form_ptr;
        best_ptm_ptr = ptm_ptr;
      }
    }
  }

  if (best_match_score == 0 || best_form_ptr == nullptr) {
    return nullptr;
  }
  else {
    double best_mass_shift = adjust_prec_mass - best_form_ptr->getMass();
    return onePtmLocalize(best_form_ptr, extend_ms_ptr_vec, best_mass_shift, 
                          best_match_score, best_ptm_ptr);
  }
}


int LocalProcessor::compOnePtmScr(ProteoformPtr base_form_ptr, 
                                  const ExtendMsPtrVec & extend_ms_ptr_vec,
                                  PtmPtr ptm_ptr) {

  int len = base_form_ptr->getLen();
  // score table with size len  
  std::vector<int> s_table;
  //LOG_DEBUG("ptm mass " << ptm_ptr->getMonoMass());
  local_util::compOnePtmSTable(s_table, base_form_ptr, extend_ms_ptr_vec, ptm_ptr, mng_ptr_);

  int max_score = 0;
  for (int i = 0; i < len; i++) {
    //LOG_DEBUG(i << " " << s_table[i]);
    if (modifiable(base_form_ptr, i, ptm_ptr) && s_table[i] > max_score) {
      max_score = s_table[i];
    }
  }
  return max_score;
}

ProteoformPtr LocalProcessor::onePtmLocalize(ProteoformPtr base_form_ptr, 
                                             const ExtendMsPtrVec & extend_ms_ptr_vec,
                                             double shift_mass, int match_score,
                                             PtmPtr ptm_ptr) {

  int len = base_form_ptr->getLen();
  // score table with size len  
  std::vector<int> s_table;
  local_util::compOnePtmSTable(s_table, base_form_ptr, extend_ms_ptr_vec, ptm_ptr, mng_ptr_);

  std::vector<double> scr_vec;
  for (int i = 0; i < len; i++) {
    if (modifiable(base_form_ptr, i, ptm_ptr)) {
      double score = std::pow(mng_ptr_->prob_ratio_, s_table[i]); 
      scr_vec.push_back(score);
    }
    else {
      scr_vec.push_back(0);
    }
  }
  local_util::normalize(scr_vec);

  ProteoformPtr local_form_ptr = local_proteoform::createProteoformPtr(base_form_ptr, match_score, 
                                                                       scr_vec, ptm_ptr, mng_ptr_); 
  return local_form_ptr;
}

bool massShiftOverlap(ProteoformPtr form_ptr_1, ProteoformPtr form_ptr_2) {
  MassShiftPtrVec shift_ptr_vec_1 = form_ptr_1->getMassShiftPtrVec(AlterType::UNEXPECTED);
  std::sort(shift_ptr_vec_1.begin(), shift_ptr_vec_1.end(), MassShift::cmpPosInc);
  MassShiftPtrVec shift_ptr_vec_2 = form_ptr_2->getMassShiftPtrVec(AlterType::VARIABLE);
  std::sort(shift_ptr_vec_2.begin(), shift_ptr_vec_2.end(), MassShift::cmpPosInc);
  if (shift_ptr_vec_1.size() != 2 || shift_ptr_vec_2.size() != 2) {
    LOG_ERROR("Mass shift number is incorrect!");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < 2; i++) {
    MassShiftPtr shift_1 = shift_ptr_vec_1[i];
    MassShiftPtr shift_2 = shift_ptr_vec_2[i];
    int left_1 = form_ptr_1->getStartPos() + shift_1->getLeftBpPos();
    int right_1 = form_ptr_1->getStartPos() + shift_1->getRightBpPos();
    int left_2 = form_ptr_2->getStartPos() + shift_2->getLeftBpPos();
    int right_2 = form_ptr_2->getStartPos() + shift_2->getRightBpPos();
    if (right_1 <= left_2 || right_2 <= left_1) {
      return false;
    }
  }
  return true;
}


// similar to processOneKnownPtm, we might get a nullptr from this function
ProteoformPtr LocalProcessor::processTwoKnownPtms(PrsmPtr prsm_ptr) {
  ProteoformPtr ori_form_ptr = prsm_ptr->getProteoformPtr();
  ProteoformPtrVec cand_form_vec = local_proteoform::getAllCandidateForms(ori_form_ptr, mng_ptr_);

  ExtendMsPtrVec extend_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  double adjust_prec_mass = prsm_ptr->getAdjustedPrecMass();
  double err_tole = mng_ptr_->peak_tole_ptr_->compStrictErrorTole(adjust_prec_mass);
  
  // get best form and best ptm pair
  int best_match_score = 0;
  ProteoformPtr best_form_ptr = nullptr;
  PtmPair best_ptm_pair;

  for (size_t i = 0; i < cand_form_vec.size(); i++) {
    ProteoformPtr cand_form_ptr = cand_form_vec[i];
    // Check if the unexpected mass shift can be explained by a PTM
    // sum of all expected shift values, in most cases, there is only one unexpected shift
    double unexp_shift_mass = adjust_prec_mass - cand_form_ptr->getMass(); 
    //LOG_DEBUG("Get PTM by Mass");
    // Get candidate Ptm pairs with similar mass shifts
    PtmPairVec ptm_pair_vec = local_util::getPtmPairVecByMass(unexp_shift_mass, err_tole, ptm_pair_vec_);

    for (size_t k = 0; k < ptm_pair_vec.size(); k++) {
      int score = compTwoPtmScr(cand_form_ptr, extend_ms_ptr_vec,
                                ptm_pair_vec[k].first, ptm_pair_vec[k].second);
      if (score > best_match_score) {
        best_match_score = score;
        best_form_ptr = cand_form_ptr;
        best_ptm_pair = ptm_pair_vec[k];
      }
    }
  }

  LOG_DEBUG("best_match_score " << best_match_score);
  if (best_form_ptr != nullptr) {
    LOG_DEBUG(best_form_ptr->getProteoformMatchSeq());
  }
  if (best_ptm_pair.first != nullptr) {
    LOG_DEBUG(best_ptm_pair.first->getName());
  }
  if (best_ptm_pair.second != nullptr) {
    LOG_DEBUG(best_ptm_pair.second->getName());
  }


  int ori_match_score = prsm_ptr->getMatchFragNum();
  if (best_match_score <= 0 
      || best_match_score < ori_match_score - mng_ptr_->DESC_MATCH_LIMIT_
      || best_match_score < ori_match_score * mng_ptr_->desc_ratio_) {
    return nullptr;
  }
  else {
    ProteoformPtr result_ptr = twoPtmLocalize(best_form_ptr, extend_ms_ptr_vec, best_match_score, 
                                              best_ptm_pair.first, best_ptm_pair.second);
    // check if the PTM localization results are overlapping with the original
    // mass shifts
    
    //result_ptr can be null when the proteoform with ptm yields low score than filtering threshold
    //in that case, running massShiftOverlap causes seg fault
    if (result_ptr == nullptr) {
      return nullptr;
    }
    if (massShiftOverlap(ori_form_ptr, result_ptr)) {
      return result_ptr;
    }
    else {
      return nullptr;
    }
  }
}

int LocalProcessor::compTwoPtmScr(ProteoformPtr form_ptr, const ExtendMsPtrVec & extend_ms_ptr_vec,
                                  PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2) {
  int len = form_ptr->getLen();
  if (len <= 0) {
    return 0;
  }

  // matching table with size 3 * (len + 1)
  std::vector<std::vector<int>> m_table;
  local_util::compTwoPtmMTable(m_table, form_ptr, extend_ms_ptr_vec, ptm_ptr_1, ptm_ptr_2, mng_ptr_);

  // fill D(f, g)
  int d_table[3][len + 1];
  memset(d_table, 0, sizeof(int) * 3 * (len + 1));
  d_table[0][0] = 0;
  d_table[1][0] = std::numeric_limits<int>::min();
  d_table[2][0] = std::numeric_limits<int>::min();

  // layer 0
  for (int i = 1; i <= len; i++) {
    d_table[0][i] = d_table[0][i-1] + m_table[0][i];
    //LOG_DEBUG("row 1 " << i << " " << d_table[0][i]);
  }

  // layer 1
  for (int i = 1; i <= len; i++) {
    d_table[1][i] = d_table[1][i-1] + m_table[1][i];
    if (modifiable(form_ptr, i - 1, ptm_ptr_1) && d_table[0][i-1] > d_table[1][i-1]) {
      //LOG_DEBUG("MODIFY" << i );
      d_table[1][i] = d_table[0][i-1] + m_table[1][i];
    }
    //LOG_DEBUG("row 2 " << i << " " << d_table[1][i]);
  }

  // layer 2
  for (int i = 1; i <= len; i++) {
    d_table[2][i] = d_table[2][i-1] + m_table[2][i];
    if (modifiable(form_ptr, i - 1, ptm_ptr_2) && d_table[1][i-1] > d_table[2][i-1]) {
      d_table[2][i] = d_table[1][i-1] + m_table[2][i];
    } 
  }

  return d_table[2][len]; 
}

  
ProteoformPtr LocalProcessor::twoPtmLocalize(ProteoformPtr form_ptr, const ExtendMsPtrVec &extend_ms_ptr_vec,
                                             int match_score, PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2) {
  int len = form_ptr->getLen();
  // score table with size 3 * (len + 1)
  std::vector<std::vector<int>> m_table;
  local_util::compTwoPtmMTable(m_table, form_ptr, extend_ms_ptr_vec, ptm_ptr_1, ptm_ptr_2, mng_ptr_);
  std::vector<int> n_term_scores = local_util::compPrefScore(m_table[0]); 
  std::vector<int> c_term_scores = local_util::compSuffScore(m_table[2]); 

  int mid_scores[len+1][len + 1];
  memset(mid_scores, 0, sizeof(int) * (len+1) * (len + 1));
  for (int bgn = 0; bgn < len + 1; bgn ++) {
    int total_score = 0;
    for (int end = bgn ; end < len + 1; end++) {
      total_score = total_score + m_table[1][end]; 
      mid_scores[bgn][end] = total_score;
      //LOG_DEBUG("mid score " << bgn << " " << end << " " << total_score);
    }
  }
  std::vector<int> ptm_1_pos_list;
  std::vector<int> ptm_2_pos_list;
  for (int i = 0; i < len; i++) {
    if (modifiable(form_ptr, i, ptm_ptr_1)) {
      ptm_1_pos_list.push_back(i);
      LOG_DEBUG("ptm 1 pos " << i);
    }
    if (modifiable(form_ptr, i, ptm_ptr_2)) {
      ptm_2_pos_list.push_back(i);
    }
  }

  double b_table[len][len];
  for (int i = 0; i < len; i++) {
    for (int j = 0; j < len; j++) {
      b_table[i][j] = 0.0;
    }
  }
  double total_prob = 0.0;
  for (size_t i = 0; i < ptm_1_pos_list.size(); i++) {
    int pos_1 = ptm_1_pos_list[i];
    for (size_t j = 0; j < ptm_2_pos_list.size(); j++) {
      int pos_2 = ptm_2_pos_list[j];
      if (pos_2 > pos_1) {
        int score = n_term_scores[pos_1] + c_term_scores[pos_2 + 1] + mid_scores[pos_1+1][pos_2];
        b_table[pos_1][pos_2] = std::pow(mng_ptr_->prob_ratio_, score); 
        total_prob = total_prob + b_table[pos_1][pos_2];
        //LOG_DEBUG("pos " << pos_1 << " " << pos_2 << " " << score << " " << b_table[pos_1][pos_2] << " " << total_prob);
      }
    }
  }

  if (total_prob == 0.0) {
    LOG_ERROR("No valid modification sites!");
    exit(EXIT_FAILURE);
  }

  // normalize 
  for (size_t i = 0; i < ptm_1_pos_list.size(); i++) {
    int pos_1 = ptm_1_pos_list[i];
    for (size_t j = 0; j < ptm_2_pos_list.size(); j++) {
      int pos_2 = ptm_2_pos_list[j];
      if (pos_2 > pos_1) {
        b_table[pos_1][pos_2] = b_table[pos_1][pos_2]/total_prob;
      }
    }
  }
  
  // find breaking points
  std::vector<double> break_scores;
  double s = 0; 
  // break point is the right size of the position
  for (int pos = 0; pos < len - 1; pos++) {
    // row to add 
    int row = pos;
    for (int j = row + 1; j < len; j++) {
      s = s +  b_table[row][j];
    }
    // column to remove
    int col = pos;
    for (int i = 0; i < col; i++) {
      s = s - b_table[i][col];
    }
    break_scores.push_back(s);
    //LOG_DEBUG("break score " << pos << " " << s);
  }
  int best_break_pos = std::distance(break_scores.begin(), std::max_element(break_scores.begin(), break_scores.end()));
  LOG_DEBUG("break pos " << best_break_pos);

  // get two score vectors
  std::vector<double> scr_vec_1 (len, 0);
  for (size_t i = 0; i < ptm_1_pos_list.size(); i++) {
    int pos_1 = ptm_1_pos_list[i];
    if (pos_1 > best_break_pos) {
      break;
    }
    double pos_scr = 0;
    for (int j = ptm_2_pos_list.size() - 1; j >= 0;  j--) {
      int pos_2 = ptm_2_pos_list[j];
      //LOG_DEBUG("pos 1 " << pos_1 << " pos 2 " << pos_2 << " " << b_table[pos_1][pos_2] << " " << pos_scr);
      if (pos_2 > pos_1 && pos_2 > best_break_pos) {
        pos_scr = pos_scr + b_table[pos_1][pos_2];
      }
      else {
        break;
      }
    }
    scr_vec_1[pos_1] = pos_scr;
    //LOG_DEBUG("scr vec one  " << pos_1 << " " << pos_scr);
  }
  std::vector<double> scr_vec_2 (len, 0);
  for (size_t j = 0; j < ptm_2_pos_list.size(); j++) {
    int pos_2 = ptm_2_pos_list[j];
    if (pos_2 <= best_break_pos) {
      continue; 
    }
    double pos_scr = 0;
    for (size_t i = 0; i < ptm_1_pos_list.size(); i++) {
      int pos_1 = ptm_1_pos_list[i];
      if (pos_1 < pos_2 && pos_1 <= best_break_pos) {
        pos_scr = pos_scr + b_table[pos_1][pos_2];
      }
      else {
        break;
      }
    }
    scr_vec_2[pos_2] = pos_scr;
    //LOG_DEBUG("scr vec two  " << pos_2 << " " << pos_scr);
  }

  ProteoformPtr local_form_ptr = local_proteoform::createProteoformPtr(form_ptr, match_score, best_break_pos, 
                                                                       scr_vec_1, scr_vec_2, ptm_ptr_1, ptm_ptr_2, mng_ptr_); 

  return local_form_ptr;
}

}  // namespace toppic
