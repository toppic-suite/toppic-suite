//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#ifndef TOPPIC_SEARCH_ONE_PTM_SEARCH_DIAGONAL_HPP_
#define TOPPIC_SEARCH_ONE_PTM_SEARCH_DIAGONAL_HPP_

#include <memory>
#include <vector>

#include "spec/theo_peak.hpp"
#include "spec/extend_ms.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "search/oneptmsearch/diagonal_header.hpp"

namespace toppic {

template <class T>
class Diagonal {
 public:
  Diagonal() {}

  Diagonal(DiagonalHeaderPtr header_ptr) {
    header_ptr_ = header_ptr; 
  }

  // need init pair_ptr_list after create
  Diagonal(DiagonalHeaderPtr header_ptr, std::vector<T> pair_ptr_list) {
    header_ptr_ = header_ptr;
    pair_ptr_list_ = pair_ptr_list;
  }

  size_t size() {return pair_ptr_list_.size(); }

  DiagonalHeaderPtr getHeader() {return header_ptr_;}

  const std::vector<T>& getDiagPair() {return pair_ptr_list_;}

  T getDiagPair(int i) {return pair_ptr_list_[i];}

 private:
  DiagonalHeaderPtr header_ptr_;

  std::vector<T> pair_ptr_list_;
};

double refinePrecursorAndHeaderShift(ProteoformPtr proteo_ptr,
                                     const ExtendMsPtrVec &ms_three_ptr_vec,
                                     DiagonalHeaderPtrVec &header_ptrs,
                                     double ppo, double min_mass,
                                     double refine_prec_step_width);

DiagonalHeaderPtrVec refineHeadersBgnEnd(ProteoformPtr proteo_ptr,
                                         const ExtendMsPtrVec &ms_three_ptr_vec,
                                         const DiagonalHeaderPtrVec& heade_ptrs,
                                         double min_mass);

DiagonalHeaderPtrVec2D refineHeadersBgnEnd(ProteoformPtr proteo_ptr,
                                           const ExtendMsPtrVec &ms_three_ptr_vec,
                                           const DiagonalHeaderPtrVec2D& header_ptrs_2d,
                                           const DiagonalHeaderPtrVec& header_ptrs_1d,
                                           double min_mass);

int getNewBgn(const PeakIonPairPtrVec& pair_ptrs);

int getNewEnd(const PeakIonPairPtrVec& pair_ptrs);

} /* namespace toppic */

#endif /* DIAGONAL_HPP_ */
