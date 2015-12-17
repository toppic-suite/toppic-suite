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
#include "ptmsearch/basic_diag_pair.hpp"
#include "ptmsearch/ps_align.hpp"
//#include "ptmsearch/diagonal.hpp"

namespace prot {

class PtmSlowMatch {
 public:
  PtmSlowMatch(ProteoformPtr proteo_ptr,
      SpectrumSetPtr spectrum_set_ptr,
      CompShiftLowMemPtr comp_shift_ptr,
      PtmMngPtr mng_ptr);

  ProteoformPtr getProteoform(){return proteo_ptr_;};
  void compute(AlignTypePtr type_ptr, PrsmPtrVec &prsm_ptrs);

 private:
  PtmMngPtr mng_ptr_;
  ProteoformPtr proteo_ptr_;
  double prec_mono_mass_;
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  PrmMsPtrVec ms_six_ptr_vec_;
  ExtendMsPtrVec ms_three_ptr_vec_;
  PSAlignPtr ps_align_ptr_;

  void initPsAlign(CompShiftLowMemPtr comp_shift_ptr);

  DiagonalHeaderPtrVec getNTermShiftHeaders(
      const std::vector<double> &best_shifts, double prec_mono_mass, 
      ProteoformPtr proteo_ptr, PtmMngPtr mng_ptr);

  PrsmPtr geneResult(int shift_num);
};

typedef std::shared_ptr<PtmSlowMatch> PtmSlowMatchPtr;
typedef std::vector<PtmSlowMatchPtr> PtmSlowMatchPtrVec;

} /* namespace prot */

#endif /* PTM_SLOW_MATCH_HPP_ */
