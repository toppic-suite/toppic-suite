#ifndef PROT_ONE_PTM_SLOW_MATCH_HPP_
#define PROT_ONE_PTM_SLOW_MATCH_HPP_

#include <memory>
#include <vector>

#include "base/proteoform.hpp"
#include "spec/prm_peak.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "oneptmsearch/diagonal_header.hpp"
#include "oneptmsearch/one_ptm_search_mng.hpp"
#include "oneptmsearch/comp_shift_low_mem.hpp"
#include "oneptmsearch/one_ptm_search_mng.hpp"
#include "oneptmsearch/basic_diag_pair.hpp"
#include "oneptmsearch/ps_align.hpp"

namespace prot {

class OnePtmSlowMatch {
 public:
  OnePtmSlowMatch(ProteoformPtr proteo_ptr,
                  SpectrumSetPtr spectrum_set_ptr,
                  CompShiftLowMemPtr comp_shift_ptr,
                  OnePtmSearchMngPtr mng_ptr);

  ProteoformPtr getProteoform(){return proteo_ptr_;};
  void compute(AlignTypePtr type_ptr, PrsmPtrVec &prsm_ptrs);

 private:
  OnePtmSearchMngPtr mng_ptr_;
  ProteoformPtr proteo_ptr_;
  double prec_mono_mass_;
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  PrmMsPtrVec ms_six_ptr_vec_;
  ExtendMsPtrVec ms_three_ptr_vec_;
  PSAlignPtr ps_align_ptr_;

  void initPsAlign(CompShiftLowMemPtr comp_shift_ptr);

  DiagonalHeaderPtrVec getNTermShiftHeaders(
      const std::vector<double> &best_shifts, double prec_mono_mass, 
      ProteoformPtr proteo_ptr, OnePtmSearchMngPtr mng_ptr);

  PrsmPtr geneResult(int shift_num);
};

typedef std::shared_ptr<OnePtmSlowMatch> OnePtmSlowMatchPtr;
typedef std::vector<OnePtmSlowMatchPtr> OnePtmSlowMatchPtrVec;

} /* namespace prot */

#endif 
