#ifndef PROT_ZERO_PTM_SEARCH_HPP_
#define PROT_ZERO_PTM_SEARCH_HPP_

#include <array>

#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmsearch/zero_ptm_mng.hpp"

namespace prot {


//void zeroPtmSearchAll(SpectrumSetPtr spectrum_set_ptr,  std::array<SimplePrSMPtrVec> &prsms);

void zeroPtmSearch(SpectrumSetPtr spec_set_ptr, int type,
                   ProteoformPtrVec &form_ptr_vec, int report_num,
                   SimplePrSMPtrVec &prsms);

void zeroPtmSearchProcess(ZeroPtmMngPtr);

} /* namespace_prot */

#endif
