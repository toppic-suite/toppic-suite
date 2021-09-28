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
#include "ms/factory/extend_ms_util.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"

#include "prsm/prsm.hpp"
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

  //logger::log_level=2; 
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
          LOG_DEBUG("Start localization");
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

  /*
  // Check if the mass shift can be explained by two common PTMs
  ProteoformPtr two_known_proteoform = processTwoKnownPtm(prsm);

  if (two_known_proteoform != nullptr) {
  double new_num_match_ion = local_util::compMatchFragNum(two_known_proteoform,
  prsm->getRefineMsPtrVec(),
  mng_ptr_->min_mass_);
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

/*
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
   */


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
  double best_mass_shift = 0;
  ProteoformPtr best_form_ptr = nullptr;
  PtmPtr best_ptm_ptr = nullptr;

  for (size_t i = 0; i < cand_form_vec.size(); i++) {
    ProteoformPtr cand_form_ptr = cand_form_vec[i];
    double unexp_shift_mass = adjust_prec_mass - cand_form_ptr->getMass(); 
    // Get candidate Ptms with similar mass shifts
    PtmPtrVec match_ptm_ptr_vec = local_util::getPtmPtrVecByMass(unexp_shift_mass, err_tole, ptm_ptr_vec_);

    // if there is a match, find the best ptm and its best similarity score
    for (size_t j = 0; j < match_ptm_ptr_vec.size(); j++) {
      PtmPtr ptm_ptr = match_ptm_ptr_vec[j];
      // compute similarity score for each possible site of the PTM
      int score = compOnePtmScr(cand_form_ptr, extend_ms_ptr_vec, ptm_ptr); 
      if (score > best_match_score) {
        best_match_score = score;
        best_mass_shift = unexp_shift_mass;
        best_form_ptr = cand_form_ptr;
        best_ptm_ptr = ptm_ptr;
      }
    }
  }

  if (best_match_score == 0 || best_form_ptr == nullptr) {
    return nullptr;
  }
  else {
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
  local_util::compOnePtmSTable(s_table, base_form_ptr, extend_ms_ptr_vec, ptm_ptr, mng_ptr_);

  int max_score = 0;
  for (int i = 0; i < len; i++) {
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

  ProteoformPtr local_form_ptr = local_proteoform::createProteoformPtr(base_form_ptr, shift_mass, match_score, 
                                                                       scr_vec, ptm_ptr, mng_ptr_); 
  return local_form_ptr;
}


// similar to processOneKnownPtm, we might get a nullptr from this function
ProteoformPtr LocalProcessor::processTwoKnownPtm(PrsmPtr prsm_ptr) {
  ProteoformPtr ori_form_ptr = prsm_ptr->getProteoformPtr();
  ProteoformPtrVec cand_form_vec = local_proteoform::getAllCandidateForms(ori_form_ptr, mng_ptr_);

  ExtendMsPtrVec extend_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  double adjust_prec_mass = prsm_ptr->getAdjustedPrecMass();
  double err_tole = mng_ptr_->peak_tole_ptr_->compStrictErrorTole(adjust_prec_mass);
  
  // get best form and best ptm pair
  int best_match_score = 0;
  double best_mass_shift = 0.0;
  ProteoformPtr best_form_ptr = nullptr;
  PtmPair best_ptm_pair;

  for (size_t i = 0; i < cand_form_vec.size(); i++) {
    ProteoformPtr cand_form_ptr = cand_form_vec[i];
    // Check if the unexpected mass shift can be explained by a PTM
    // sum of all expected shift values, in most cases, there is only one unexpected shift
    double unexp_shift_mass = adjust_prec_mass - cand_form_ptr->getMass(); 
    LOG_DEBUG("Get PTM by Mass");
    // Get candidate Ptm pairs with similar mass shifts
    PtmPairVec ptm_pair_vec = local_util::getPtmPairVecByMass(unexp_shift_mass, err_tole, ptm_pair_vec_);

    for (size_t k = 0; k < ptm_pair_vec.size(); k++) {
      int score = 0;
      // TO DO
      // int score = compTwoPtmScr(cand_form_ptr, extend_ms_ptr_vec,
      //                          ptm_pair_vec[k].first, ptm_pair_vec[k].second);
      if (score > best_match_score) {
        best_match_score = score;
        best_mass_shift = unexp_shift_mass;
        best_form_ptr = cand_form_ptr;
        best_ptm_pair = ptm_pair_vec[k];
      }
    }
  }

  int ori_match_score = prsm_ptr->getMatchPeakNum();

  if (best_match_score == 0 
      || best_match_score < ori_match_score - mng_ptr_->DESC_MATCH_LIMIT_
      || best_match_score < ori_match_score * mng_ptr_->desc_ratio_) {
    return nullptr;
  }

  // TO DO localization
 
  //ProteoformPtr result_ptr = twoPtmLocalize(best_form_ptr, extend_ms_ptr_vec, prec_mass, err_tole);
  //if (result_ptr != nullptr) {
  //  return result_ptr;
  //}
  return nullptr;
}

void compTwoPtmMTable(std::vector<std::vector<int>> &s_table, ProteoformPtr form_ptr, 
                      const ExtendMsPtrVec & extend_ms_ptr_vec,
                      PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2, LocalMngPtr mng_ptr) {

  int g = form_ptr->getLen();
  for (int i = 0; i < 3; i++) {
    std::vector<int> row(g + 1, 0);
    s_table.push_back(row);
  }

  BpSpecPtr bp_spec_ptr = form_ptr->getBpSpecPtr();
  std::vector<double> prm_masses = bp_spec_ptr->getPrmMasses();
  std::vector<double> srm_masses = bp_spec_ptr->getSrmMasses();
  // sort srm in the decreasing order
  std::sort(srm_masses.begin(), srm_masses.end(), std::greater<double>());
  PeakTolerancePtr tole_ptr = mng_ptr->peak_tole_ptr_;

  int mass_1 = ptm_ptr_1->getMonoMass();
  int mass_2 = ptm_ptr_2->getMonoMass();

  for (size_t i = 0; i < extend_ms_ptr_vec.size(); i++) {
    std::vector<double> ms_masses = extend_ms_util::getExtendMassVec(extend_ms_ptr_vec[i]);
    // updated S table using prm masses 
    double n_shift = extend_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getNShift();
    local_util::compMTable(prm_masses, n_shift, ms_masses, tole_ptr, s_table[0]);  
    local_util::updateSTable(prm_masses, n_shift + mass_1, ms_masses, tole_ptr, s_table[1]);  
    local_util::updateSTable(prm_masses, n_shift + mass_1 + mass_2, ms_masses, tole_ptr, s_table[2]);  

    // updated S table using srm masses 
    double c_shift = extend_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getCShift();
    local_util::updateSTable(srm_masses, c_shift + mass_1 + mass_2, ms_masses, tole_ptr, s_table[0]);  
    local_util::updateSTable(srm_masses, c_shift + mass_2, ms_masses, tole_ptr, s_table[1]);  
    local_util::updateSTable(srm_masses, c_shift, ms_masses, tole_ptr, s_table[2]);  
  }
}

double LocalProcessor::compTwoPtmScr(ProteoformPtr form_ptr, const ExtendMsPtrVec & extend_ms_ptr_vec,
                                     PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2) {
  int g = form_ptr->getLen();
  if (g <= 0) {
    return 0;
  }

  // matching table with size 3 * (g + 1)
  std::vector<std::vector<int>> m_table;
  compTwoPtmMTable(m_table, form_ptr, extend_ms_ptr_vec, ptm_ptr_1, ptm_ptr_2, mng_ptr_);

  // fill D(f, g)
  int d_table[3][g + 1];
  memset(d_table, 0, sizeof(int) * 3 * (g + 1));
  d_table[0][0] = 1;

  // layer 0
  for (int i = 1; i <= g; i++) {
    d_table[0][i] = d_table[0][i-1] + m_table[0][i];
  }

  // layer 1
  for (int i = 1; i <= g; i++) {
    d_table[1][i] = d_table[1][i-1] + m_table[1][i];
    if (modifiable(form_ptr, i - 1, ptm_ptr_1) && d_table[0][i-1] > d_table[1][i-1]) {
      d_table[1][i] = d_table[0][i-1] + m_table[1][i];
    }
  }

  // layer 2
  for (int i = 1; i <= g; i++) {
    d_table[2][i] = d_table[2][i-1] + m_table[2][i];
    if (modifiable(form_ptr, i - 1, ptm_ptr_2) && d_table[1][i-1] > d_table[2][i-1]) {
      d_table[2][i] = d_table[1][i-1] + m_table[2][i];
    } 
  }

  return d_table[2][g]; 
}

  
/*
ProteoformPtr LocalProcessor::twoPtmLocalize(ProteoformPtr form_ptr, const ExtendMsPtrVec &extend_ms_ptr_vec,
                                             PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2) {
  // score table with size 3 * (g + 1)
  std::vector<std::vector<int>> s_table;
  compTwoPtmSTable(s_table, form_ptr, extend_ms_ptr_vec, ptm_ptr_1, ptm_ptr_2, mng_ptr_);
  int len = form_ptr->getLen();
  std::vector<int> n_term_scores (len+1, 0);
  int zero_shift_total_score = 0;
  for (int i = 0; i < len+1; i++) {
    zero_shift_total_score = zero_shift_total_score + s_table[0][i];
    n_term_scores[i] = zero_shift_total_score;
  }

  std::vector<int> c_term_scores (len+1, 0);
  int two_shift_total_score = 0;
  for (int i = len; i >=0; i--) {
    two_shift_total_score = two_shift_total_score + s_table[2][i];
    c_term_scores[i] = two_shift_total_score;
  }

  std::vector<int> ptm_1_pos_list;
  std::vector<int> ptm_2_pos_list;
  for (int i = 0; i < len; i++) {
    if (modifiable(form_ptr, i, ptm_ptr_1)) {
      ptm_1_pos_list.push_back(i);
    }
    if (modifiable(form_ptr, i, ptm_ptr_2)) {
      ptm_2_pos_list.push_back(i);
    }
  }

  int a_table[len][len];
  std::fill((int *)a_table, (int *)a_table + (len*len), 0);
  double b_table[len][len];
  double zero = 0.0;
  std::fill((double *)b_table, (double *)b_table + (len*len), zero);

  for (size_t i = 0; i < ptm_1_pos_list.size(); i++) {
    int pos_1 = ptm_1_pos_list[i];
    for (size_t j = 0; j < ptm_2_pos_list.size(); j++) {
      int pos_2 = ptm_2_pos_list[j];
      if (pos_2 > pos_1) {
        int score = n_term_scores[pos_1] + c_term_scores[pos_2 + 1] + mid_scores[pos_1+1, pos_2];
        a_table[pos_1][pos_2] = score;
        b_table[pos_1][pos_2] =
      }
    }
  }


  return nullptr;
}
*/




}  // namespace toppic

