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

#ifndef TOPPIC_SEARCH_DIAG_DIAG_PAIR_HPP_
#define TOPPIC_SEARCH_DIAG_DIAG_PAIR_HPP_

#include <memory>
#include <vector>
#include "search/diag/pair.hpp"

namespace toppic {

class DiagPair;
typedef std::shared_ptr<DiagPair> DiagPairPtr;
typedef std::vector<DiagPairPtr>  DiagPairPtrVec;

class Diagonal;
typedef std::shared_ptr<Diagonal>   DiagonalPtr;
typedef std::weak_ptr<Diagonal>   DiagonalWeakPtr;

class DiagPair:public Pair {
  public:
  DiagPair(int x, int y, double score, int diag_order, double diff):
    Pair(x, y),
    score_(score),
    diag_order_(diag_order),
    diff_(diff) {}

  int getDiagOrder() {return diag_order_;}

  const DiagonalWeakPtr geneDiagonalPtr() {return diagonal_ptr_;}

  void setDiagonalPtr(DiagonalPtr diagonal_ptr) {diagonal_ptr_ = diagonal_ptr;}

  double getDiff() {return diff_;}

  double getScore() {return score_;}

 protected:
  double score_;

  int diag_order_;

  double diff_;

  DiagonalWeakPtr diagonal_ptr_;
};

} /* namespace toppic */

#endif /* DIAG_PAIR_HPP_ */
