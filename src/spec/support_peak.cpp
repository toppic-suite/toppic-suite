
#include <spec/support_peak.hpp>

namespace prot {

SupportPeak::SupportPeak(DeconvPeakPtr deconv_peak_ptr, double offset,
                         double score, SPTypePtr peak_type_ptr){
  deconv_peak_ptr_ = deconv_peak_ptr;
  offset_=offset;
  score_=score;
  peak_type_ptr_ = peak_type_ptr;
}

} /* namespace prot */
