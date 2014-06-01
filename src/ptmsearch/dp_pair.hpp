/*
 * dp_pair.hpp
 *
 *  Created on: Jan 8, 2014
 *      Author: xunlikun
 */

#ifndef PROT_DP_PAIR_HPP_
#define PROT_DP_PAIR_HPP_

#include "ptmsearch/diagonal_header.hpp"
#include "prsm/pair.hpp"

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
         int order,int n_shift,DiagonalHeaderPtr header);

  const DPPairPtr& getDiagPrev() const {
    return diag_prev_;
  }

  void setDiagPrev(const DPPairPtr& diagPrev) {
    diag_prev_ = diagPrev;
  }

  double getDiff() const {
    return diff_;
  }

  const DiagonalHeaderPtr& getDiagonalHeader() const {
    return header_;
  }

  int getDiagOrder() const {
    return order_;
  }

  double getPairScore() const {
    return pair_score_;
  }

  double getScr(int s) {
    return scores_[s];
  }

  DPPairPtr getPre(int s){
    return prevs_[s];
  }


  int getType(int s){
    return types_[s];
  }

  bool isAssisting(){
    if(pair_score_==0.0){
      return true;
    }
    return false;
  }

  void updateTable(int s,double score,int path_type,DPPairPtr prev_pair);

 private:
  DiagonalHeaderPtr header_;
  double diff_;
  double pair_score_;
  int order_;
  DPPairPtr diag_prev_;
  DPPairPtrVec prevs_;
  std::vector<double> scores_;
  std::vector<int> types_;
};

} /* namespace prot */

#endif /* DP_PAIR_HPP_ */
