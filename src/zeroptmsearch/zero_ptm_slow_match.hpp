#ifndef ZERO_PTM_SLOW_MATCH_HPP_
#define ZERO_PTM_SLOW_MATCH_HPP_

#include "spec/deconv_peak.hpp"
#include "spec/theo_peak.hpp"
#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class ZeroPtmSlowMatch {
 public:
  ZeroPtmSlowMatch(DeconvMsPtr deconv_ms_ptr, ZpFastMatchPtr fast_match_ptr,
                   ZeroPtmMngPtr mng_ptr);

  double getScore() {return score_;}

  PrsmPtr geneResult();

 private:
  ZeroPtmMngPtr mng_ptr_;
  ZpFastMatchPtr fast_match_ptr_;
  DeconvMsPtr deconv_ms_ptr_;
  ProteoformPtr proteoform_ptr_;

  double refine_prec_mass_;
  ExtendMsPtr refine_ms_ptr_;

  double score_ = 0;
  double recal_ = 0;

  void compScore (ExtendMsPtr refine_ms_ptr, const TheoPeakPtrVec &theo_peak_ptrs, 
                  double ppo);
  bool isValid (double recal, double ppo);
};

typedef std::shared_ptr<ZeroPtmSlowMatch> ZpSlowMatchPtr;
typedef std::vector<ZpSlowMatchPtr> ZpSlowMatchPtrVec;

inline bool compareZeroPtmSlowMatchDown(const ZpSlowMatchPtr &a, 
                                        const ZpSlowMatchPtr &b) {
  return a->getScore() > b->getScore();
}

ZpSlowMatchPtrVec zeroPtmSlowFilter(DeconvMsPtr deconv_ms_ptr,
                                    const ZpFastMatchPtrVec &fast_match_ptrs,
                                    ZeroPtmMngPtr mng_ptr); 

}

#endif
