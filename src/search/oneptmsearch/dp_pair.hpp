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

#ifndef TOPPIC_SEARCH_ONE_PTM_SEARCH_DP_PAIR_HPP_
#define TOPPIC_SEARCH_ONE_PTM_SEARCH_DP_PAIR_HPP_

#include "search/diag/pair.hpp"
#include "search/diag/diag_header.hpp"

namespace toppic {

class DpPair;
typedef std::shared_ptr<DpPair>  DpPairPtr;
typedef std::vector<DpPairPtr> DpPairPtrVec;

class DpPair : public Pair{
 public:
  DpPair(int x,int y,double pair_score,double diff,
         int order,int n_shift, DiagHeaderPtr header_ptr);

  DpPairPtr getDiagPrevPairPtr() {return diag_prev_pair_ptr_;}

  void setDiagPrevPairPtr(DpPairPtr diag_prev_pair_ptr) {
    diag_prev_pair_ptr_ = diag_prev_pair_ptr;}

  double getDiff() {return diff_;}

  DiagHeaderPtr getDiagHeader() {return header_ptr_;}

  int getDiagOrder() {return order_;}

  double getPairScore() {return pair_score_;}

  double getScore(int s) {return scores_[s];}

  DpPairPtr getPrevPairPtr(int s){return prev_pair_ptrs_[s];}

  int getType(int s){return types_[s];}

  bool isAssisting(){ return (pair_score_==0.0); }

  void updateTable(int s,double score,int path_type,DpPairPtr prev_pair);

 private:
  DiagHeaderPtr header_ptr_;
  double diff_;
  double pair_score_;
  int order_;
  DpPairPtr diag_prev_pair_ptr_;
  DpPairPtrVec prev_pair_ptrs_;
  std::vector<double> scores_;
  std::vector<int> types_;
};

} /* namespace toppic */

#endif /* DP_PAIR_HPP_ */
