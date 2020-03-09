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
#include <iostream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"

#include "seq/proteoform_util.hpp"

#include "filter/mng/topindex_file_name.hpp"

#include "filter/massmatch/filter_protein.hpp"
#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "filter/zeroptm/mass_zero_ptm_filter.hpp"

namespace toppic {

MassZeroPtmFilter::MassZeroPtmFilter(const ProteoformPtrVec &proteo_ptrs,
                                     ZeroPtmFilterMngPtr mng_ptr, std::string block_str) {
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  
  std::string index_dir = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx";

	TopIndexFileNamePtr file_name_ptr = std::make_shared<TopIndexFileName>();
  std::string parameters = file_name_ptr->geneFileName(prsm_para_ptr);
  std::string suffix = parameters + block_str;//last part of file name

  //check if all index files for this ptm is present. if not, generate index files again.

 bool index_files_exist = true;

  for (size_t t = 0; t < file_name_ptr->zero_ptm_file_vec_.size(); t++){
    std::string file_name = file_name_ptr->zero_ptm_file_vec_[t] + suffix;
    if (!file_util::exists(index_dir + file_util::getFileSeparator() + file_name)){
      index_files_exist = false;//if any of the index files for this ptm is missing
      break; 
    }
  }

  if (index_files_exist){
    std::cout << "Loading index files                            " << std::endl;

    term_index_ptr_ = std::make_shared<MassMatch>();
    diag_index_ptr_ = std::make_shared<MassMatch>();
    rev_term_index_ptr_ = std::make_shared<MassMatch>();
    rev_diag_index_ptr_ = std::make_shared<MassMatch>();

    term_index_ptr_->deserializeMassMatch(file_name_ptr->zero_ptm_file_vec_[0] + suffix, index_dir);
    diag_index_ptr_->deserializeMassMatch(file_name_ptr->zero_ptm_file_vec_[1] + suffix, index_dir);
    rev_term_index_ptr_->deserializeMassMatch(file_name_ptr->zero_ptm_file_vec_[2] + suffix, index_dir);
    rev_diag_index_ptr_->deserializeMassMatch(file_name_ptr->zero_ptm_file_vec_[3] + suffix, index_dir);

  }
  
  else{
    LOG_DEBUG("get shifts");
    std::vector<std::vector<double> > shift_2d
        = proteoform_util::getNTermShift2D(proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    std::vector<std::vector<double> > n_term_acet_2d
        = proteoform_util::getNTermAcet2D(proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    LOG_DEBUG("get shifts complete");
    // N-terminal indexes
    term_index_ptr_ = MassMatchFactory::getPrmTermMassMatchPtr(proteo_ptrs, shift_2d,
                                                              mng_ptr->max_proteoform_mass_,
                                                              mng_ptr->filter_scale_);
    // Prm indexes
    diag_index_ptr_ = MassMatchFactory::getPrmDiagMassMatchPtr(proteo_ptrs,
                                                              mng_ptr->max_proteoform_mass_,
                                                              mng_ptr->filter_scale_);
    LOG_DEBUG("diag index");
    std::vector<std::vector<double> > rev_shift_2d;
    std::vector<double> shift_1d(1, 0);
    for (size_t i = 0; i < proteo_ptrs.size(); i++) {
      rev_shift_2d.push_back(shift_1d);
    }
    // C-terminal indexes
    rev_term_index_ptr_ = MassMatchFactory::getSrmTermMassMatchPtr(proteo_ptrs, rev_shift_2d,
                                                                  n_term_acet_2d,
                                                                  mng_ptr->max_proteoform_mass_,
                                                                  mng_ptr->filter_scale_);

    // To generate SRM indexes, n terminal acetylation shifts are added into the SRM list. 
    rev_diag_index_ptr_ = MassMatchFactory::getSrmDiagMassMatchPtr(proteo_ptrs, n_term_acet_2d,
                                                                  mng_ptr->max_proteoform_mass_,
                                                                  mng_ptr->filter_scale_);
 }
}
void MassZeroPtmFilter::computeBestMatch(const ExtendMsPtrVec &ms_ptr_vec) {
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  bool pref = true;
  std::vector<std::pair<int, int> > pref_mass_errors
      = extend_ms::getExtendIntMassErrorList(ms_ptr_vec, pref, mng_ptr_->filter_scale_);
  pref = false;
  std::vector<std::pair<int, int> > suff_mass_errors
      = extend_ms::getExtendIntMassErrorList(ms_ptr_vec, pref, mng_ptr_->filter_scale_);
  std::pair<int, int> prec_minus_water_mass_error
      = ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWaterError(tole_ptr->getPpo(),
                                                                        mng_ptr_->filter_scale_);

  int term_row_num = term_index_ptr_->getRowNum();
  std::vector<short> term_scores(term_row_num, 0);
  term_index_ptr_->compMatchScores(pref_mass_errors, prec_minus_water_mass_error, term_scores);
  /*
     for (int i = 0; i < term_row_num; i++) {
     LOG_DEBUG("row " << i << " score "<< term_scores[i]);
     }
     */
  int diag_row_num = diag_index_ptr_->getRowNum();
  std::vector<short> diag_scores(diag_row_num, 0);
  diag_index_ptr_->compMatchScores(pref_mass_errors, prec_minus_water_mass_error, diag_scores);

  int rev_term_row_num = rev_term_index_ptr_->getRowNum();
  std::vector<short> rev_term_scores(rev_term_row_num, 0);
  rev_term_index_ptr_->compMatchScores(suff_mass_errors, prec_minus_water_mass_error, rev_term_scores);
  /*
     for (int i = 0; i < term_row_num; i++) {
     LOG_DEBUG("rev row " << i << " score "<< rev_term_scores[i]);
     }
     */
  int rev_diag_row_num = rev_diag_index_ptr_->getRowNum();
  std::vector<short> rev_diag_scores(rev_diag_row_num, 0);
  rev_diag_index_ptr_->compMatchScores(suff_mass_errors, prec_minus_water_mass_error, rev_diag_scores);

  int threshold = MassMatch::getPrecursorMatchScore() * 2 + 4;
  bool add_shifts = false;
  int shift_num = 0;
  FilterProteinPtrVec comp_prots
      = mass_match_util::findTopProteins(term_scores, rev_term_scores, 
                                         term_index_ptr_, rev_term_index_ptr_,
                                         threshold, mng_ptr_->comp_num_, add_shifts, shift_num);
  comp_match_ptrs_.clear();
  int group_spec_num = ms_ptr_vec.size();
  for (size_t i = 0; i < comp_prots.size(); i++) {
    int id = comp_prots[i]->getProteinId();
    int score = comp_prots[i]->getScore();
    comp_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id], score));
  }

  FilterProteinPtrVec pref_prots
      = mass_match_util::findTopProteins(term_scores, rev_diag_scores, 
                                         term_index_ptr_, rev_diag_index_ptr_,
                                         threshold, mng_ptr_->pref_suff_num_, add_shifts, shift_num);
  pref_match_ptrs_.clear();
  for (size_t i = 0; i < pref_prots.size(); i++) {
    int id = pref_prots[i]->getProteinId();
    int score = pref_prots[i]->getScore();
    pref_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id], score));
  }

  FilterProteinPtrVec suff_prots
      = mass_match_util::findTopProteins(diag_scores, rev_term_scores, 
                                         diag_index_ptr_, rev_term_index_ptr_,
                                         threshold, mng_ptr_->pref_suff_num_, add_shifts, shift_num);
  suff_match_ptrs_.clear();
  for (size_t i = 0; i < suff_prots.size(); i++) {
    int id = suff_prots[i]->getProteinId();
    int score = suff_prots[i]->getScore();
    suff_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id], score));
  }

  FilterProteinPtrVec internal_prots
      = mass_match_util::findTopProteins(diag_scores, rev_diag_scores, diag_index_ptr_, rev_diag_index_ptr_,
                                         threshold, mng_ptr_->inte_num_, add_shifts, shift_num);
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_prots.size(); i++) {
    int id = internal_prots[i]->getProteinId();
    int score = internal_prots[i]->getScore();
    internal_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                                group_spec_num,
                                                                proteo_ptrs_[id], score));
  }
}




} /* namespace toppic */
