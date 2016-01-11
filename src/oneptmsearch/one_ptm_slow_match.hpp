#ifndef PROT_ONE_PTM_SLOW_MATCH_HPP_
#define PROT_ONE_PTM_SLOW_MATCH_HPP_

#include <memory>
#include <vector>

#include "base/proteoform.hpp"
#include "spec/prm_peak.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "prsm/simple_prsm.hpp"
#include "oneptmsearch/diagonal_header.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"
//#include "ptmsearch/comp_shift_low_mem.hpp"
#include "oneptmsearch/basic_diag_pair.hpp"
#include "oneptmsearch/ps_align.hpp"

namespace prot {

class OnePtmSlowMatch {
 public:
  OnePtmSlowMatch(ProteoformPtr proteo_ptr,
                  SpectrumSetPtr spectrum_set_ptr,
                  SimplePrsmPtr simple_prsm_ptr,
                  AlignTypePtr align_type_ptr,
                  PtmSearchMngPtr mng_ptr);

  ProteoformPtr getProteoform(){return proteo_ptr_;};

  void init();

  PrsmPtr compute(AlignTypePtr align_type_ptr, int shift_num);

 private:
  PtmSearchMngPtr mng_ptr_;
  ProteoformPtr proteo_ptr_;
  double prec_mono_mass_;
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  PrmMsPtrVec ms_six_ptr_vec_;
  ExtendMsPtrVec ms_three_ptr_vec_;
  AlignTypePtr align_type_ptr_;
  SimplePrsmPtr simple_prsm_ptr_;
  PSAlignPtr ps_align_ptr_;

  void addPrefixDiagonals(DiagonalHeaderPtrVec &n_extend_header_ptrs);

  void addSuffixDiagonals(DiagonalHeaderPtrVec &c_extend_header_ptrs);

  void addComplementDiagonals(DiagonalHeaderPtrVec &n_extend_header_ptrs,
                              DiagonalHeaderPtrVec &c_extend_header_ptrs);

  DiagonalHeaderPtrVec geneOnePtmNTermShiftHeaders();
};

typedef std::shared_ptr<OnePtmSlowMatch> OnePtmSlowMatchPtr;
typedef std::vector<OnePtmSlowMatchPtr> OnePtmSlowMatchPtrVec;

} /* namespace prot */

#endif 
