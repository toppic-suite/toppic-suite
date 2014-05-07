/*
 * ptm_slow_match.hpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#ifndef PROT_PTM_SLOW_MATCH_HPP_
#define PROT_PTM_SLOW_MATCH_HPP_

#include <memory>
#include <vector>
#include "base/proteoform.hpp"
#include "spec/prm_peak.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/diagonal_header.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/diagonal.hpp"
#include "ptmsearch/basic_diag_pair.hpp"
#include "ptmsearch/ps_align.hpp"

namespace prot {

class PtmSlowMatch {
 public:
  PtmSlowMatch(
      ProteoformPtr seq,
      SpectrumSetPtr spectrum_set,
      CompShiftLowMemPtr comp_shift,
      PtmMngPtr mng);
  ProteoformPtr getSeq(){return seq_;};
  double getScr(int shiftnum,int type);
  PrsmPtr geneResult(int shift_num);
  void compute(SemiAlignTypePtr type, PrsmPtrVec &prsms);

 private:
  PtmMngPtr mng_;
  ProteoformPtr seq_;
  DeconvMsPtr deconv_ms_;
  PrmMsPtr ms_six_;
  ExtendMsPtr ms_three_;
  PSAlignPtr ps_align_;

  std::vector<double> scores_;

  DiagonalHeaderPtrVec2D headers_;

  DiagonalHeaderPtrVec2D prefix_headers_;
  DiagonalHeaderPtrVec2D suffix_headers_;
  DiagonalHeaderPtrVec2D internal_headers_;
  std::vector<double> complete_deltas_;
  std::vector<double> prefix_deltas_;
  std::vector<double> suffix_deltas_;
  std::vector<double> internal_deltas_;

  void init(CompShiftLowMemPtr comp_shift);

  DiagonalHeaderPtrVec getNTermShiftList(
      std::vector<double> best_shift,
      PrmMsPtr ms_six,
      ProteoformPtr seq,
      PtmMngPtr mng);

  bool found(
      double shift,
      DiagonalHeaderPtrVec headerlist,
      double trunc_error_tolerance);
};

typedef std::shared_ptr<PtmSlowMatch> PtmSlowMatchPtr;
typedef std::vector<PtmSlowMatchPtr> PtmSlowMatchPtrVec;

} /* namespace prot */

#endif /* PTM_SLOW_MATCH_HPP_ */
