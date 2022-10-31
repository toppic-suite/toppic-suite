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

#include <boost/thread/mutex.hpp>

#include "common/util/file_util.hpp"
#include "ms/factory/prm_ms_util.hpp"
#include "prsm/simple_prsm_util.hpp"

#include "filter/massmatch/prot_candidate.hpp"
#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "filter/diag/diag_filter.hpp"

namespace toppic {

// serialization mutex.
boost::mutex diag_filter_mutex;

DiagFilter::DiagFilter(const ProteoformPtrVec &proteo_ptrs,
                       DiagFilterMngPtr mng_ptr, std::string block_str) {
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  
  std::string parameters = mng_ptr->getIndexFilePara();
  std::string suffix = parameters + block_str;
	
  std::string index_dir = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx";

  bool index_files_exist = true;

  for (size_t t = 0; t < mng_ptr->multi_ptm_file_vec_.size(); t++){
    std::string file_name = mng_ptr->multi_ptm_file_vec_[t] + suffix;
    if (!file_util::exists(index_dir + file_util::getFileSeparator() + file_name)){
      //if any of the index files for this ptm is missing
      index_files_exist = false;
      break; 
    }
  }

  if (index_files_exist){ 
    index_ptr_ = std::make_shared<MassMatch>();
    std::string file_name = mng_ptr->multi_ptm_file_vec_[0] + suffix;

    {
      boost::unique_lock<boost::mutex> lock(diag_filter_mutex);
      index_ptr_->deserializeMassMatch(file_name, index_dir);
    }
  }
  else{
    index_ptr_ = mass_match_factory::getPrmDiagMassMatchPtr(proteo_ptrs,
                                                            mng_ptr->max_proteoform_mass_,
                                                            mng_ptr->filter_scale_);
  }
}

SimplePrsmPtrVec DiagFilter::getBestMatch(const PrmMsPtrVec &ms_ptr_vec) {
  SimplePrsmPtrVec match_ptrs = compute(ms_ptr_vec);
  SimplePrsmPtrVec unique_match_ptrs = simple_prsm_util::getUniqueMatches(match_ptrs);
  std::sort(unique_match_ptrs.begin(), unique_match_ptrs.end(), SimplePrsm::cmpScoreDec);
  size_t num = mng_ptr_->filter_result_num_;
  if (num > unique_match_ptrs.size()) {
    num = unique_match_ptrs.size();
  }
  SimplePrsmPtrVec result_ptrs;
  for (size_t i = 0; i < num; i++) {
    SimplePrsmPtr match_ptr = unique_match_ptrs[i];
    if (match_ptr->getScore() > 0.0) {
      result_ptrs.push_back(match_ptr);
    } else {
      break;
    }
  }
  return result_ptrs;
}

SimplePrsmPtrVec DiagFilter::compute(const PrmMsPtrVec &ms_ptr_vec) {
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int, int> > mass_errors
      = prm_ms_util::getIntMassErrorList(ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_, true, false);
  SimplePrsmPtrVec match_ptrs;
  int row_num = index_ptr_->getRowNum();
  std::vector<int> proteo_row_begins = index_ptr_->getProteoRowBegins();
  std::vector<int> proteo_row_ends = index_ptr_->getProteoRowEnds();
  int threshold = 4;
  std::vector<short> scores(row_num, 0);
  std::vector<short> max_scores(row_num, 0);
  for (size_t i = 0; i < mass_errors.size(); i++) {
    std::fill(scores.begin(), scores.end(), 0);
    index_ptr_->compScores(mass_errors, i, -mass_errors[i].first, scores);
    for (int j = 0; j < row_num; j++) {
      if (scores[j] > max_scores[j]) {
        max_scores[j] = scores[j];
      }
    }
  }
 
  ProtCandidatePtrVec results
    = mass_match_util::findTopProteins(max_scores, proteo_row_begins, proteo_row_ends, threshold,
                                       mng_ptr_->filter_result_num_);
  for (size_t j = 0; j < results.size(); j++) {
    int id = results[j]->getProteinId();
    int score = results[j]->getScore();
    match_ptrs.push_back(std::make_shared<SimplePrsm>(ms_ptr_vec[0]->getMsHeaderPtr(),
                                                      ms_ptr_vec.size(),
                                                      proteo_ptrs_[id], score));
  }
  return match_ptrs;
}

} /* namespace toppic */
