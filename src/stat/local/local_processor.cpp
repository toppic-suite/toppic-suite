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

  return prsm;
}

// we will get a nullptr if the mass shift can't be explained by a variable ptm
ProteoformPtr LocalProcessor::processOneKnownPtm(PrsmPtr prsm_ptr) {
  //get canidate forms
  ProteoformPtr ori_form_ptr = prsm_ptr->getProteoformPtr();
  ProteoformPtrVec cand_form_vec = local_proteoform::getAllCandidateForms(ori_form_ptr, mng_ptr_);

  ExtendMsPtrVec extend_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  double adjust_prec_mass = prsm_ptr->getAdjustedPrecMass();
  double err_tole = mng_ptr_->peak_tole_ptr_->compStrictErrorTole(adjust_prec_mass);
  
  int best_match_score = -1;
  ProteoformPtr best_form_ptr = nullptr;

  for (size_t i = 0; i < cand_form_vec.size(); i++) {
    ProteoformPtr cand_form_ptr = cand_form_vec[i];
    LOG_DEBUG("Start one ptm localization");
    LocalResultPtr local_result_ptr = onePtmLocalize(cand_form_ptr, extend_ms_ptr_vec, 
                                                     adjust_prec_mass, err_tole); 
    LOG_DEBUG("End one ptm localization");
    if (local_result_ptr != nullptr && local_result_ptr->match_score_ > best_match_score) {
      best_match_score = local_result_ptr->match_score_;
      best_form_ptr = local_result_ptr->form_ptr_;
    }
  }
  return best_form_ptr;
}

// similar to processOneKnownPtm, we might get a nullptr from this function
ProteoformPtr LocalProcessor::processTwoKnownPtm(PrsmPtr prsm_ptr) {
  ProteoformPtr ori_form_ptr = prsm_ptr->getProteoformPtr();
  ProteoformPtrVec cand_form_vec = local_proteoform::getAllCandidateForms(ori_form_ptr, mng_ptr_);

  ExtendMsPtrVec extend_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  double adjust_prec_mass = prsm_ptr->getAdjustedPrecMass();
  double err_tole = mng_ptr_->peak_tole_ptr_->compStrictErrorTole(adjust_prec_mass);
  
  int best_match_score = 0;
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
      int score = dpTwoPtmScr(cand_form_ptr, extend_ms_ptr_vec,
                              ptm_pair_vec[k].first, ptm_pair_vec[k].second);
      if (score > best_match_score) {
        best_match_score = score;
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

  return nullptr;
}

LocalResultPtr LocalProcessor::onePtmLocalize(ProteoformPtr form_ptr, const ExtendMsPtrVec & extend_ms_ptr_vec, 
                                              double prec_mass, double err_tole) {

  // Check if the unexpected mass shift can be explained by a PTM
  // sum of all expected shift values, in most cases, there is only one unexpected shift
  double unexp_shift_mass = prec_mass - form_ptr->getMass(); 

  LOG_DEBUG("Get PTM by Mass");
  // Get candidate Ptms with similar mass shifts
  PtmPtrVec match_ptm_ptr_vec = local_util::getPtmPtrVecByMass(unexp_shift_mass, err_tole, ptm_ptr_vec_);

  // if there is a match, find the best ptm and its best similarity score
  if (match_ptm_ptr_vec.size() != 0) {
    LOG_DEBUG("computer ptm score");
    // compute similarity score for each possible site of the PTM
    LocalResultPtr local_result_ptr = compOnePtmScr(form_ptr, extend_ms_ptr_vec, unexp_shift_mass, match_ptm_ptr_vec);
    return local_result_ptr;
  }
  return nullptr;
}

/*
LocalResultPtr LocalProcessor::twoPtmLocalize(ProteoformPtr form_ptr, const ExtendMsPtrVec & extend_ms_ptr_vec, 
                                              double prec_mass, double err_tole) {


  // if there is a match, find the best ptm and its best similarity score
  if (match_ptm_pair_vec.size() != 0) {
    LOG_DEBUG("computer ptm score");
    // compute similarity score for each possible site of the PTM
    // LocalResultPtr local_result_ptr = compTwoPtmScr(form_ptr, extend_ms_ptr_vec, unexp_shift_mass, 
    //                                                 match_ptm_pair_vec, ori_match_num);
    return local_result_ptr;
  }
  return nullptr;
}
*/

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
        unexp_shift_ptr->setMassShift(ptm_ptr_vec[i]->getMonoMass());
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
    ProteoformPtr local_form_ptr = local_proteoform::createProteoformPtr(base_form_ptr, shift_mass, 
                                                                         result_ptr, mng_ptr_);
    if (local_form_ptr != nullptr) {
      result_ptr->form_ptr_ = local_form_ptr;
      return result_ptr;
    }
  }

  return nullptr;
}

double LocalProcessor::dpTwoPtmScr(ProteoformPtr form_ptr, const ExtendMsPtrVec & extend_ms_ptr_vec,
                                   PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2) {
  int g = form_ptr->getLen();
  if (g <= 0) {
    return 0;
  }

  int mass_1 = ptm_ptr_1->getMonoMass();
  int mass_2 = ptm_ptr_2->getMonoMass();

  // score table with size 3 * g
  std::vector<std::vector<int>> s_table;
  for (int i = 0; i < 3; i++) {
    std::vector<int> row(g + 1, 0);
    s_table.push_back(row);
  }

  BpSpecPtr bp_spec_ptr = form_ptr->getBpSpecPtr();
  std::vector<double> prm_masses = bp_spec_ptr->getPrmMasses();
  std::vector<double> srm_masses = bp_spec_ptr->getSrmMasses();
  // sort srm in the decreasing order
  std::sort(srm_masses.begin(), srm_masses.end(), std::greater<double>());
  PeakTolerancePtr tole_ptr = mng_ptr_->peak_tole_ptr_;

  for (size_t i = 0; i < extend_ms_ptr_vec.size(); i++) {
    std::vector<double> ms_masses = extend_ms_util::getExtendMassVec(extend_ms_ptr_vec[i]);
    // updated S table using prm masses 
    double n_shift = extend_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getNShift();
    local_util::updateSTable(prm_masses, n_shift, ms_masses, tole_ptr, s_table[0]);  
    local_util::updateSTable(prm_masses, n_shift + mass_1, ms_masses, tole_ptr, s_table[1]);  
    local_util::updateSTable(prm_masses, n_shift + mass_1 + mass_2, ms_masses, tole_ptr, s_table[2]);  

    // updated S table using srm masses 
    double c_shift = extend_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getCShift();
    local_util::updateSTable(srm_masses, c_shift + mass_1 + mass_2, ms_masses, tole_ptr, s_table[0]);  
    local_util::updateSTable(srm_masses, c_shift + mass_2, ms_masses, tole_ptr, s_table[1]);  
    local_util::updateSTable(srm_masses, c_shift, ms_masses, tole_ptr, s_table[2]);  
  }

  // fill D(f, g)
  int d_table[3][g + 1];
  memset(d_table, 0, sizeof(int) * 3 * (g + 1));
  d_table[0][0] = 1;

  // layer 0
  for (int i = 1; i <= g; i++) {
    d_table[0][i] = d_table[0][i-1] + s_table[0][i];
  }

  // layer 1
  for (int i = 1; i <= g; i++) {
    d_table[1][i] = d_table[1][i-1] + s_table[1][i];
    if (modifiable(form_ptr, i - 1, ptm_ptr_1) && d_table[0][i-1] > d_table[1][i-1]) {
      d_table[1][i] = d_table[0][i-1] + s_table[1][i];
    }
  }

  // layer 2
  for (int i = 1; i <= g; i++) {
    d_table[2][i] = d_table[2][i-1] + s_table[2][i];
    if (modifiable(form_ptr, i - 1, ptm_ptr_2) && d_table[1][i-1] > d_table[2][i-1]) {
      d_table[2][i] = d_table[1][i-1] + s_table[2][i];
    } 
  }

  return d_table[2][g]; 
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

