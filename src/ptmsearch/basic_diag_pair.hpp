#ifndef PROT_BASIC_DIAG_PAIR_HPP_
#define PROT_BASIC_DIAG_PAIR_HPP_

#include <memory>
#include <vector>
#include "base/algorithm.hpp"
#include "spec/prm_peak.hpp"
#include "ptmsearch/pair.hpp"
#include "ptmsearch/diagonal.hpp"
#include "ptmsearch/diagonal_header.hpp"

namespace prot {

class BasicDiagPair;
typedef std::shared_ptr<BasicDiagPair> BasicDiagPairPtr;
typedef std::vector<BasicDiagPairPtr> BasicDiagPairPtrVec;
typedef Diagonal<BasicDiagPairPtr> BasicDiagonal;
typedef std::shared_ptr<BasicDiagonal> BasicDiagonalPtr;
typedef std::weak_ptr<BasicDiagonal> BasicDiagonalWeakPtr;
typedef std::vector<BasicDiagonalPtr> BasicDiagonalPtrVec;

class BasicDiagPair:public Pair {
 public:
  BasicDiagPair(int x,int y,double score,int diag_order, double diff);

  int getDiagOrder() {return diag_order_;}

  const BasicDiagonalWeakPtr geneDiagonalPtr() {return diagonal_ptr_;}

  void setDiagonalPtr(BasicDiagonalPtr diagonal_ptr) {diagonal_ptr_ = diagonal_ptr;}

  double getDiff() {return diff_;}

  double getScore() {return score_;}

 protected:
  int diag_order_;
  double diff_;
  BasicDiagonalWeakPtr diagonal_ptr_;
  double score_;
};


BasicDiagonalPtrVec geneDiagonals(const DiagonalHeaderPtrVec& header_ptr_vec,
                                 const PrmPeakPtrVec &prm_peaks, 
                                 int group_spec_num, ProteoformPtr proteo_ptr);
/*
BasicDiagPairPtrVec compDiagPair(const PrmMsPtrVec &ms_ptr_vec, 
                                 const std::vector<double>& seq_masses,
                                 DiagonalHeaderPtr header_ptr);

BasicDiagonalPtr  getDiagonalPtr(int cnt, DiagonalHeaderPtr header_ptr, 
                                 const PrmMsPtrVec &ms_six_ptr_vec,
                                 ProteoformPtr proteo_ptr, PtmMngPtr mng_ptr);
                                 */
} /* namespace prot */

#endif /* BASIC_DIAG_PAIR_HPP_ */
