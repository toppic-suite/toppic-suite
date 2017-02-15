// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_DP_PAIR_HPP_
#define PROT_DP_PAIR_HPP_

#include "oneptmsearch/pair.hpp"
#include "oneptmsearch/diagonal_header.hpp"

namespace prot {

#define PATH_TYPE_NULL -1
#define PATH_TYPE_DIAGONAL 0
#define PATH_TYPE_SHIFT 1
#define PATH_TYPE_TRUNC 2

class DPPair;
typedef std::shared_ptr<DPPair>  DPPairPtr;
typedef std::vector<DPPairPtr> DPPairPtrVec;

class DPPair : public Pair{
 public:
  DPPair(int x,int y,double pair_score,double diff,
         int order,int n_shift, DiagonalHeaderPtr header_ptr);

  DPPairPtr getDiagPrevPairPtr() {return diag_prev_pair_ptr_;}

  void setDiagPrevPairPtr(DPPairPtr diag_prev_pair_ptr) {
    diag_prev_pair_ptr_ = diag_prev_pair_ptr;}

  double getDiff() {return diff_;}

  DiagonalHeaderPtr getDiagonalHeader() {return header_ptr_;}

  int getDiagOrder() {return order_;}

  double getPairScore() {return pair_score_;}

  double getScore(int s) {return scores_[s];}

  DPPairPtr getPrevPairPtr(int s){return prev_pair_ptrs_[s];}

  int getType(int s){return types_[s];}

  bool isAssisting(){ return (pair_score_==0.0); }

  void updateTable(int s,double score,int path_type,DPPairPtr prev_pair);

 private:
  DiagonalHeaderPtr header_ptr_;
  double diff_;
  double pair_score_;
  int order_;
  DPPairPtr diag_prev_pair_ptr_;
  DPPairPtrVec prev_pair_ptrs_;
  std::vector<double> scores_;
  std::vector<int> types_;
};

} /* namespace prot */

#endif /* DP_PAIR_HPP_ */
