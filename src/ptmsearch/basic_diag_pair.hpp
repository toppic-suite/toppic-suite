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

  int getBaseType() {return base_type_;}

  int getDiagOrder() {return diag_order_;}

  const BasicDiagPairDiagPtr& getDiagonal() {return diagonal_;}

  void setDiagonal(const BasicDiagPairDiagPtr& diagonal) {diagonal_ = diagonal;}

  double getDiff() {return diff_;}

  double getScore() {return score_;}

 protected:
  int diag_order_;
  double diff_;
  BasicDiagPairDiagPtr diagonal_;
  double score_;
  int base_type_;

};

BasicDiagPairPtrVec compDiagPair(PrmMsPtr ms_ptr, const std::vector<double>& seq_masses,
                                 DiagonalHeaderPtr header_ptr);

BasicDiagPairDiagPtrVec getDiagonals(const DiagonalHeaderPtrVec& header_ptrs,
                                     PrmMsPtr ms_six_ptr, ProteoformPtr proteo_ptr,
                                     PtmMngPtr mng_ptr);

BasicDiagPairDiagPtr getDiagonal(int cnt,DiagonalHeaderPtr header_ptr, PrmMsPtr ms_six_ptr,
                                 ProteoformPtr proteo_ptr, PtmMngPtr mng_ptr);
} /* namespace prot */

#endif /* BASIC_DIAG_PAIR_HPP_ */
