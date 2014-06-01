/*
 * basic_diag_pair.hpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

#ifndef PROT_BASIC_DIAG_PAIR_HPP_
#define PROT_BASIC_DIAG_PAIR_HPP_

#include <memory>
#include <vector>
#include "base/algorithm.hpp"
#include "spec/prm_peak.hpp"
#include "prsm/pair.hpp"
#include "ptmsearch/diagonal.hpp"
#include "ptmsearch/diagonal_header.hpp"

namespace prot {

class BasicDiagPair;
typedef std::shared_ptr<BasicDiagPair> BasicDiagPairPtr;
typedef std::vector<BasicDiagPairPtr> BasicDiagPairPtrVec;
typedef std::shared_ptr<Diagonal<BasicDiagPairPtr>> BasicDiagPairDiagPtr;
typedef std::vector<BasicDiagPairDiagPtr> BasicDiagPairDiagPtrVec;

class BasicDiagPair:public Pair {
 public:
  BasicDiagPair(int x,int y,double score,int diag_order,
                double diff,int prm_peak_type);

  int getBaseType() const {
    return base_type_;
  }

  int getDiagOrder() const {
    return diag_order_;
  }

  const BasicDiagPairDiagPtr& getDiagonal() const {
    return diagonal_;
  }

  void setDiagonal(const BasicDiagPairDiagPtr& diagonal) {
    diagonal_ = diagonal;
  }

  double getDiff() const {
    return diff_;
  }

  double getScore() const {
    return score_;
  }

 protected:
  int diag_order_;
  double diff_;
  BasicDiagPairDiagPtr diagonal_;
  double score_;
  int base_type_;

};

bool contains(BasicDiagPairPtrVec pairs,int y);

BasicDiagPairPtrVec compDiagPair(PrmMsPtr sp,std::vector<double> seq_masses,
                                 DiagonalHeaderPtr header);

BasicDiagPairDiagPtrVec getDiagonals(DiagonalHeaderPtrVec headers,
                                     PrmMsPtr ms_six,ProteoformPtr seq,
                                     PtmMngPtr mng);

BasicDiagPairDiagPtr getDiagonal(int cnt,DiagonalHeaderPtr header,PrmMsPtr ms_six,
                                 ProteoformPtr seq,PtmMngPtr mng);
} /* namespace prot */

#endif /* BASIC_DIAG_PAIR_HPP_ */
