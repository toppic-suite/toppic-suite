#ifndef PROT_ONE_PTM_SLOW_MATCH_HPP_
#define PROT_ONE_PTM_SLOW_MATCH_HPP_

#include <memory>
#include <vector>

#include "base/proteoform.hpp"
#include "spec/prm_peak.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/diagonal_header.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"
#include "ptmsearch/basic_diag_pair.hpp"
#include "ptmsearch/ps_align.hpp"

#include "oneptmsearch/one_ptm_search_mng.hpp"

namespace prot {

class OnePtmSlowMatch {
 public:
  OnePtmSlowMatch(ProteoformPtr proteo_ptr,
                  SpectrumSetPtr spectrum_set_ptr,
                  SemiAlignTypePtr align_type_ptr,
                  OnePtmSearchMngPtr mng_ptr);

  ProteoformPtr getProteoform(){return proteo_ptr_;};
  PrsmPtr getResultPrsm();

 private:
  OnePtmSearchMngPtr mng_ptr_;
  SemiAlignTypePtr align_type_ptr_;
  ProteoformPtr proteo_ptr_;
  double prec_mono_mass_;
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  PrmMsPtrVec ms_six_ptr_vec_;
  ExtendMsPtrVec ms_three_ptr_vec_;
  PrmPeakPtrVec prm_peaks_;
  int group_spec_num_;

  PSAlignPtr ps_align_ptr_;

  void initOnePtmAlign();

  BasicDiagonalPtrVec compAlignDiagonals();
  std::vector<double> compBestShifts(CompShiftLowMemPtr comp_shift_ptr);
  void extendNHeaders(DiagonalHeaderPtrVec &n_extend_header_ptrs);
  void extendCHeaders(DiagonalHeaderPtrVec &c_extend_header_ptrs);
  BasicDiagonalPtrVec removeEmptyDiagonals(BasicDiagonalPtrVec &n_diagonal_ptrs,
                                           BasicDiagonalPtrVec &c_diagonal_ptrs);
  BasicDiagonalPtrVec getDiagonals(DiagonalHeaderPtrVec &header_ptrs);
};

typedef std::shared_ptr<OnePtmSlowMatch> OnePtmSlowMatchPtr;
typedef std::vector<OnePtmSlowMatchPtr> OnePtmSlowMatchPtrVec;

} /* namespace prot */

#endif /* ONE_PTM_SLOW_MATCH_HPP_ */
