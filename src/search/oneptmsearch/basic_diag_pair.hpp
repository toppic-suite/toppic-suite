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


#ifndef TOPPIC_SEARCH_ONE_PTM_SEARCH_BASIC_DIAG_PAIR_HPP_
#define TOPPIC_SEARCH_ONE_PTM_SEARCH_BASIC_DIAG_PAIR_HPP_

#include <memory>
#include <vector>
#include "ms/spec/prm_peak.hpp"
#include "search/oneptmsearch/pair.hpp"
#include "search/oneptmsearch/diagonal.hpp"
#include "search/oneptmsearch/diagonal_header.hpp"

namespace toppic {

class BasicDiagPair;
typedef std::shared_ptr<BasicDiagPair> BasicDiagPairPtr;
typedef std::vector<BasicDiagPairPtr>  BasicDiagPairPtrVec;
typedef Diagonal<BasicDiagPairPtr>     BasicDiagonal;
typedef std::shared_ptr<BasicDiagonal> BasicDiagonalPtr;
typedef std::weak_ptr<BasicDiagonal>   BasicDiagonalWeakPtr;
typedef std::vector<BasicDiagonalPtr>  BasicDiagonalPtrVec;

class BasicDiagPair:public Pair {
 public:
  BasicDiagPair(int x, int y, double score, int diag_order, double diff):
      Pair(x, y),
      score_(score),
      diag_order_(diag_order),
      diff_(diff) {}

  int getDiagOrder() {return diag_order_;}

  const BasicDiagonalWeakPtr geneDiagonalPtr() {return diagonal_ptr_;}

  void setDiagonalPtr(BasicDiagonalPtr diagonal_ptr) {diagonal_ptr_ = diagonal_ptr;}

  double getDiff() {return diff_;}

  double getScore() {return score_;}

 protected:
  double score_;

  int diag_order_;

  double diff_;

  BasicDiagonalWeakPtr diagonal_ptr_;
};

BasicDiagonalPtrVec geneDiagonals(const DiagonalHeaderPtrVec& header_ptr_vec,
                                  const PrmPeakPtrVec &prm_peaks,
                                  int group_spec_num, ProteoformPtr proteo_ptr);
} /* namespace toppic */

#endif /* BASIC_DIAG_PAIR_HPP_ */
