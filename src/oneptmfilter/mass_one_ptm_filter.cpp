#include <algorithm>
#include <iostream>
#include "base/proteoform_util.hpp"
#include "spec/prm_ms.hpp"
#include "zeroptmfilter/filter_protein.hpp"
#include "zeroptmfilter/mass_match_factory.hpp"
#include "zeroptmfilter/mass_match_util.hpp"
#include "oneptmfilter/mass_one_ptm_filter.hpp"

namespace prot {

MassOnePtmFilter::MassOnePtmFilter(const ProteoformPtrVec &proteo_ptrs,
                           OnePtmFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  double rev = false;
  std::vector<std::vector<double>> shift_2d 
      = ProteoformUtil::getNTermShift2D(proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
  term_index_ptr_ = MassMatchFactory::getMassMatchPtr(proteo_ptrs, shift_2d,
                                               mng_ptr->max_proteoform_mass_,
                                               mng_ptr->filter_scale_, rev);
  diag_index_ptr_ = MassMatchFactory::getMassMatchPtr(proteo_ptrs, 
                                                      mng_ptr->max_proteoform_mass_,
                                                      mng_ptr->filter_scale_, rev);
  rev = true;
  std::vector<std::vector<double>> rev_shift_2d;
  std::vector<double> shift_1d(1,0);
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    rev_shift_2d.push_back(shift_1d);
  }
  rev_term_index_ptr_ = MassMatchFactory::getMassMatchPtr(proteo_ptrs, rev_shift_2d,
                                                          mng_ptr->max_proteoform_mass_,
                                                          mng_ptr->filter_scale_, rev);
  rev_diag_index_ptr_ = MassMatchFactory::getMassMatchPtr(proteo_ptrs, mng_ptr->max_proteoform_mass_,
                                                          mng_ptr->filter_scale_, rev);
}

void MassOnePtmFilter::computeBestMatch(const PrmMsPtrVec &prm_ms_ptr_vec, 
                                    const PrmMsPtrVec &srm_ms_ptr_vec){
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int,int>> pref_mass_errors 
      = PrmMs::getIntMassErrorList(prm_ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_, true, false);
  std::vector<std::pair<int,int>> suff_mass_errors 
      = PrmMs::getIntMassErrorList(srm_ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_, false, true);

  int term_row_num = term_index_ptr_->getRowNum();
  std::vector<short> term_scores(term_row_num, 0);
  term_index_ptr_->compScores(pref_mass_errors, term_scores);

  int diag_row_num = diag_index_ptr_->getRowNum();
  std::vector<short> diag_scores(diag_row_num, 0);
  diag_index_ptr_->compScores(pref_mass_errors, diag_scores);

  int rev_term_row_num = rev_term_index_ptr_->getRowNum();
  std::vector<short> rev_term_scores(rev_term_row_num, 0);
  rev_term_index_ptr_->compScores(suff_mass_errors, rev_term_scores);

  int rev_diag_row_num = rev_diag_index_ptr_->getRowNum();
  std::vector<short> rev_diag_scores(rev_diag_row_num, 0);
  rev_diag_index_ptr_->compScores(suff_mass_errors, rev_diag_scores);

  int threshold = 4;
  bool add_shifts = true;
  FilterProteinPtrVec comp_prots 
      = MassMatchUtil::findTopProteins(term_scores, rev_term_scores, term_index_ptr_, rev_term_index_ptr_,
                                       threshold, mng_ptr_->comp_num_, add_shifts, mng_ptr_->shift_num_);  
  comp_match_ptrs_.clear();
  int group_spec_num = prm_ms_ptr_vec.size();
  for (size_t i = 0; i < comp_prots.size(); i++) {
    int id = comp_prots[i]->getProteinId();
    int score = comp_prots[i]->getScore();
    comp_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(prm_ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[id], score)));
  }

  FilterProteinPtrVec pref_prots 
      = MassMatchUtil::findTopProteins(term_scores, rev_diag_scores, term_index_ptr_, rev_diag_index_ptr_,
                                       threshold, mng_ptr_->pref_suff_num_, add_shifts, mng_ptr_->shift_num_);  
  pref_match_ptrs_.clear();
  for (size_t i = 0; i < pref_prots.size(); i++) {
    int id = pref_prots[i]->getProteinId();
    int score = pref_prots[i]->getScore();
    SimplePrsmPtr prsm_ptr(new SimplePrsm(prm_ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                          proteo_ptrs_[id], score));
    prsm_ptr->setCTruncShifts(pref_prots[i]->getCTermShifts());
    pref_match_ptrs_.push_back(prsm_ptr); 
  }

  FilterProteinPtrVec suff_prots 
      = MassMatchUtil::findTopProteins(diag_scores, rev_term_scores, diag_index_ptr_, rev_term_index_ptr_,
                                       threshold, mng_ptr_->pref_suff_num_, add_shifts, mng_ptr_->shift_num_);  
  suff_match_ptrs_.clear();
  for (size_t i = 0; i < suff_prots.size(); i++) {
    int id = suff_prots[i]->getProteinId();
    int score = suff_prots[i]->getScore();
    SimplePrsmPtr prsm_ptr(new SimplePrsm(prm_ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                          proteo_ptrs_[id], score));
    prsm_ptr->setNTruncShifts(suff_prots[i]->getNTermShifts());
    suff_match_ptrs_.push_back(prsm_ptr); 
  }

  FilterProteinPtrVec internal_prots 
      = MassMatchUtil::findTopProteins(diag_scores, rev_diag_scores, diag_index_ptr_, rev_diag_index_ptr_,
                                       threshold, mng_ptr_->inte_num_, add_shifts, mng_ptr_->shift_num_);  
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_prots.size(); i++) {
    int id = internal_prots[i]->getProteinId();
    int score = internal_prots[i]->getScore();
    SimplePrsmPtr prsm_ptr(new SimplePrsm(prm_ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                          proteo_ptrs_[id], score));
    prsm_ptr->setNTruncShifts(internal_prots[i]->getNTermShifts());
    prsm_ptr->setCTruncShifts(internal_prots[i]->getCTermShifts());
    internal_match_ptrs_.push_back(prsm_ptr); 
  }
}


} /* namespace prot */
