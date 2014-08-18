#ifndef PROT_PS_ALIGN_HPP_
#define PROT_PS_ALIGN_HPP_

#include "ptmsearch/dp_pair.hpp"
#include "ptmsearch/basic_diag_pair.hpp"

namespace prot {

class PSAlign {
 public:
  PSAlign();
  PSAlign(const std::vector<double> &ms_masses,
          const std::vector<double> &seq_masses,
          const BasicDiagonalPtrVec &diagonal_ptrs,
          PtmMngPtr mng_ptr);
  void compute(SemiAlignTypePtr align_type_ptr);
  void initDPPair();
  void dp(SemiAlignTypePtr align_type_ptr);
  void backtrace();
  DiagonalHeaderPtrVec backtrace(int s);

  double getAlignScr(int s){return align_scores_[s];};
  DiagonalHeaderPtrVec getResult(int s){return backtrack_diagonal_ptrs_[s];};

 protected:
  PtmMngPtr mng_;
  std::vector<double> ms_masses_;
  std::vector<double> seq_masses_;
  BasicDiagonalPtrVec diagonal_ptrs_;
  std::vector<std::vector<int>> idxes_;
  std::vector<std::vector<bool>> penalties_;

  std::vector<DPPairPtrVec> dp_2d_pair_ptrs_;
  DPPairPtr first_pair_ptr_;
  DPPairPtr last_pair_ptr_;
  DPPairPtrVec segment_bgn_pair_ptrs_;
  DPPairPtrVec segment_end_pair_ptrs_;
  DPPairPtrVec dp_pair_ptrs_;
  DiagonalHeaderPtrVec2D backtrack_diagonal_ptrs_;
  std::vector<double> align_scores_;

  void dpPrep();
  DPPairPtr getTruncPre(DPPairPtr cur_pair_ptr,int s, SemiAlignTypePtr type_ptr);
  DPPairPtr getShiftPre(DPPairPtr cur_pair_ptr,int p,int s,SemiAlignTypePtr type_ptr);
  DPPairPtr oldGetShiftPre(DPPairPtr cur_pair_ptr,int p,int s,SemiAlignTypePtr type_ptr);
};

typedef std::shared_ptr<PSAlign> PSAlignPtr;
typedef std::vector<PSAlignPtr> PSAlignPtrVec;
} /* namespace prot */

#endif /* PS_ALIGN_HPP_ */
