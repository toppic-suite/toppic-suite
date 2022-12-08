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

#include <boost/thread/mutex.hpp>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"

#include "seq/proteoform_util.hpp"
#include "seq/prot_candidate.hpp"
#include "ms/factory/prm_ms_util.hpp"

#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "filter/varptm/var_ptm_filter.hpp"

namespace toppic {

// serialization mutex.
boost::mutex var_ptm_filter_mutex;

VarPtmFilter::VarPtmFilter(const ProteoformPtrVec &proteo_ptrs,
                           VarPtmFilterMngPtr mng_ptr, 
                           std::string block_str) {
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  
  std::string index_dir = mng_ptr_->prsm_para_ptr_->getDbIndexDir(); 

  std::string parameters = mng_ptr->getIndexFilePara();
  std::string suffix = parameters + block_str;//last part of file name

  //check if all index files for this ptm is present. if not, generate index files again.
  bool index_files_exist = true;
  for (size_t t = 0; t < mng_ptr->var_ptm_file_vec_.size(); t++){
    std::string file_name = mng_ptr->var_ptm_file_vec_[t] + suffix;
    if (!file_util::exists(index_dir + file_util::getFileSeparator() + file_name)){
      //if any of the index files for this ptm is missing
      index_files_exist = false;
      break; 
    }
  }

  if (index_files_exist){
    term_index_ptr_ = std::make_shared<MassMatch>();
    diag_index_ptr_ = std::make_shared<MassMatch>();
    rev_term_index_ptr_ = std::make_shared<MassMatch>();
    rev_diag_index_ptr_ = std::make_shared<MassMatch>();
    {
      LOG_DEBUG("Loading index files");
      boost::unique_lock<boost::mutex> lock(var_ptm_filter_mutex);
      term_index_ptr_->deserializeMassMatch(mng_ptr->var_ptm_file_vec_[0] + suffix, index_dir);
      diag_index_ptr_->deserializeMassMatch(mng_ptr->var_ptm_file_vec_[1] + suffix, index_dir);
      rev_term_index_ptr_->deserializeMassMatch(mng_ptr->var_ptm_file_vec_[2] + suffix, index_dir);
      rev_diag_index_ptr_->deserializeMassMatch(mng_ptr->var_ptm_file_vec_[3] + suffix, index_dir);
    }
  }
  else{
    LOG_DEBUG("get shifts");
    std::vector<std::vector<double> > shift_2d
        = proteoform_util::getNTermShift2D(proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    std::vector<std::vector<double> > n_term_acet_2d
        = proteoform_util::getNTermAcet2D(proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
    LOG_DEBUG("get shifts complete");
    // N-terminal indexes
    term_index_ptr_ = mass_match_factory::getPrmTermMassMatchPtr(proteo_ptrs, shift_2d,
                                                                 mng_ptr->max_proteoform_mass_,
                                                                 mng_ptr->filter_scale_);
    LOG_DEBUG("diag index");
    // Prm indexes
    diag_index_ptr_ = mass_match_factory::getPrmDiagMassMatchPtr(proteo_ptrs,
                                                                 mng_ptr->max_proteoform_mass_,
                                                                 mng_ptr->filter_scale_);
    std::vector<std::vector<double> > rev_shift_2d;
    std::vector<double> shift_1d(1, 0);
    for (size_t i = 0; i < proteo_ptrs.size(); i++) {
      rev_shift_2d.push_back(shift_1d);
    }
    // C-terminal indexes
    rev_term_index_ptr_ = mass_match_factory::getSrmTermMassMatchPtr(proteo_ptrs, rev_shift_2d,
                                                                     n_term_acet_2d,
                                                                     mng_ptr->max_proteoform_mass_,
                                                                     mng_ptr->filter_scale_);

    // To generate SRM indexes, n terminal acetylation shifts are added into the SRM list. 
    rev_diag_index_ptr_ = mass_match_factory::getSrmDiagMassMatchPtr(proteo_ptrs, n_term_acet_2d,
                                                                     mng_ptr->max_proteoform_mass_,
                                                                     mng_ptr->filter_scale_);
 }
}

std::vector<std::pair<int, int>> getShiftedMassErrors(std::pair<int, int> mass_error, 
                                                      std::vector<int> shifts) {
  std::vector<std::pair<int,int>> result;
  for (size_t i = 0; i < shifts.size(); i++) {
    // skip shift 0
    if (shifts[i] == 0) {
      continue;
    }
    std::pair<int, int> m_e (mass_error.first - shifts[i], mass_error.second);
    result.push_back(m_e);
  }
  return result;
}

void VarPtmFilter::computeBestMatch(const PrmMsPtrVec &prm_ms_ptr_vec,
                                    const PrmMsPtrVec &srm_ms_ptr_vec) {
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int, int>> pref_mass_errors
      = prm_ms_util::getIntMassErrorList(prm_ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_, true, false);
  std::vector<std::pair<int, int>> suff_mass_errors
      = prm_ms_util::getIntMassErrorList(srm_ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_, false, true);
  std::pair<int, int> prec_minus_water_mass_error
    = prm_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWaterError(tole_ptr->getPpo(),
                                                                          mng_ptr_->filter_scale_);
  std::vector<std::pair<int, int> > prec_minus_water_mass_errors 
    = getShiftedMassErrors(prec_minus_water_mass_error, mng_ptr_->int_shift_list_);
  
  int term_row_num = term_index_ptr_->getRowNum();
  std::vector<short> term_scores(term_row_num, 0);
  term_index_ptr_->compMatchScores(pref_mass_errors, prec_minus_water_mass_errors, term_scores);

  int diag_row_num = diag_index_ptr_->getRowNum();
  std::vector<short> diag_scores(diag_row_num, 0);
  diag_index_ptr_->compMatchScores(pref_mass_errors, prec_minus_water_mass_errors, diag_scores);

  int rev_term_row_num = rev_term_index_ptr_->getRowNum();
  std::vector<short> rev_term_scores(rev_term_row_num, 0);
  rev_term_index_ptr_->compMatchScores(suff_mass_errors, prec_minus_water_mass_errors, rev_term_scores);

  int rev_diag_row_num = rev_diag_index_ptr_->getRowNum();
  std::vector<short> rev_diag_scores(rev_diag_row_num, 0);
  rev_diag_index_ptr_->compMatchScores(suff_mass_errors, prec_minus_water_mass_errors, rev_diag_scores);

  int group_spec_num = prm_ms_ptr_vec.size();
  double prec_minus_water_mass = prm_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWater();
  double prec_error_tole = prm_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecErrorTolerance(tole_ptr->getPpo());

  ProtCandidatePtrVec comp_prots
    = mass_match_util::findVarPtmTopProteins(term_scores, rev_term_scores, 
                                             term_index_ptr_, rev_term_index_ptr_,
                                             prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_, 
                                             mng_ptr_->comp_threshold_, 
                                             mng_ptr_->comp_num_);
  comp_match_ptrs_.clear();
  for (size_t i = 0; i < comp_prots.size(); i++) {
    int id = comp_prots[i]->getProteinId();
    comp_match_ptrs_.push_back(std::make_shared<SimplePrsm>(prm_ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id], 
                                                            comp_prots[i]));
  }

  ProtCandidatePtrVec pref_prots
    = mass_match_util::findVarPtmTopProteins(term_scores, rev_diag_scores, 
                                             term_index_ptr_, rev_diag_index_ptr_,
                                             prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             mng_ptr_->pref_suff_threshold_, 
                                             mng_ptr_->pref_suff_num_);
  pref_match_ptrs_.clear();
  for (size_t i = 0; i < pref_prots.size(); i++) {
    int id = pref_prots[i]->getProteinId();
    pref_match_ptrs_.push_back(std::make_shared<SimplePrsm>(prm_ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id],
                                                            pref_prots[i]));
  }

  ProtCandidatePtrVec suff_prots
    = mass_match_util::findVarPtmTopProteins(diag_scores, rev_term_scores, 
                                             diag_index_ptr_, rev_term_index_ptr_,
                                             prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             mng_ptr_->pref_suff_threshold_, 
                                             mng_ptr_->pref_suff_num_);
  suff_match_ptrs_.clear();
  for (size_t i = 0; i < suff_prots.size(); i++) {
    int id = suff_prots[i]->getProteinId();
    suff_match_ptrs_.push_back(std::make_shared<SimplePrsm>(prm_ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id], 
                                                            suff_prots[i]));
  }

  ProtCandidatePtrVec internal_prots
    = mass_match_util::findVarPtmTopProteins(diag_scores, rev_diag_scores, 
                                             diag_index_ptr_, rev_diag_index_ptr_,
                                             prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             mng_ptr_->internal_threshold_, 
                                             mng_ptr_->internal_num_);  
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_prots.size(); i++) {
    int id = internal_prots[i]->getProteinId();
    internal_match_ptrs_.push_back(std::make_shared<SimplePrsm>(prm_ms_ptr_vec[0]->getMsHeaderPtr(),
                                                                group_spec_num,
                                                                proteo_ptrs_[id], 
                                                                internal_prots[i]));
  }

}

/*
void VarPtmFilter::computeBestMatch(const ExtendMsPtrVec &ms_ptr_vec) {

  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  bool pref = true;
  std::vector<std::pair<int, int> > pref_mass_errors
      = extend_ms_util::getExtendIntMassErrorList(ms_ptr_vec, pref, mng_ptr_->filter_scale_);
  pref = false;
  std::vector<std::pair<int, int> > suff_mass_errors
      = extend_ms_util::getExtendIntMassErrorList(ms_ptr_vec, pref, mng_ptr_->filter_scale_);
  std::pair<int, int> prec_minus_water_mass_error
      = ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWaterError(tole_ptr->getPpo(),
                                                                        mng_ptr_->filter_scale_);
  std::vector<std::pair<int, int> > prec_minus_water_mass_errors 
    = getShiftedMassErrors(prec_minus_water_mass_error, mng_ptr_->int_shift_list_);

  int term_row_num = term_index_ptr_->getRowNum();
  std::vector<short> term_scores(term_row_num, 0);
  term_index_ptr_->compMatchScores(pref_mass_errors, prec_minus_water_mass_errors, term_scores);

  int diag_row_num = diag_index_ptr_->getRowNum();
  std::vector<short> diag_scores(diag_row_num, 0);
  diag_index_ptr_->compMatchScores(pref_mass_errors, prec_minus_water_mass_errors, diag_scores);

  int rev_term_row_num = rev_term_index_ptr_->getRowNum();
  std::vector<short> rev_term_scores(rev_term_row_num, 0);
  rev_term_index_ptr_->compMatchScores(suff_mass_errors, prec_minus_water_mass_errors, rev_term_scores);

  int rev_diag_row_num = rev_diag_index_ptr_->getRowNum();
  std::vector<short> rev_diag_scores(rev_diag_row_num, 0);
  rev_diag_index_ptr_->compMatchScores(suff_mass_errors, prec_minus_water_mass_errors, rev_diag_scores);

  int threshold = MassMatch::getPrecursorMatchScore() * 2 + 4;
  double prec_minus_water_mass = ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWater();
  double prec_error_tole = ms_ptr_vec[0]->getMsHeaderPtr()->getPrecErrorTolerance(tole_ptr->getPpo());
  int group_spec_num = ms_ptr_vec.size();
  LOG_DEBUG("mass " << ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass());
  ProtCandidatePtrVec comp_prots
    = mass_match_util::findVarPtmTopProteins(term_scores, rev_term_scores, 
                                             term_index_ptr_, rev_term_index_ptr_,
                                             prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_, 
                                             threshold, mng_ptr_->comp_num_);
  comp_match_ptrs_.clear();
  for (size_t i = 0; i < comp_prots.size(); i++) {
    int id = comp_prots[i]->getProteinId();
    comp_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id], 
                                                            comp_prots[i]));
  }

  ProtCandidatePtrVec pref_prots
    = mass_match_util::findVarPtmTopProteins(term_scores, rev_diag_scores, 
                                             term_index_ptr_, rev_diag_index_ptr_,
                                             prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             threshold, mng_ptr_->pref_suff_num_);
  pref_match_ptrs_.clear();
  for (size_t i = 0; i < pref_prots.size(); i++) {
    int id = pref_prots[i]->getProteinId();
    pref_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id],
                                                            pref_prots[i]));
  }

  ProtCandidatePtrVec suff_prots
    = mass_match_util::findVarPtmTopProteins(diag_scores, rev_term_scores, 
                                             diag_index_ptr_, rev_term_index_ptr_,
                                             prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             threshold, mng_ptr_->pref_suff_num_);
  suff_match_ptrs_.clear();
  for (size_t i = 0; i < suff_prots.size(); i++) {
    int id = suff_prots[i]->getProteinId();
    suff_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                            group_spec_num,
                                                            proteo_ptrs_[id], 
                                                            suff_prots[i]));
  }

  ProtCandidatePtrVec internal_prots
    = mass_match_util::findVarPtmTopProteins(diag_scores, rev_diag_scores, 
                                             diag_index_ptr_, rev_diag_index_ptr_,
                                             prec_minus_water_mass, prec_error_tole,
                                             mng_ptr_->shift_list_,
                                             threshold, mng_ptr_->inte_num_);  
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_prots.size(); i++) {
    int id = internal_prots[i]->getProteinId();
    internal_match_ptrs_.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                                group_spec_num,
                                                                proteo_ptrs_[id], 
                                                                internal_prots[i]));
  }
}
*/

} /* namespace toppic */
