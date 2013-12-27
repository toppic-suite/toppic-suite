#ifndef PROT_ZERO_PTM_SEARCH_HPP_
#define PROT_ZERO_PTM_SEARCH_HPP_

#include <array>

#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "zeroptmsearch/zero_ptm_mng.hpp"

namespace prot {

void zeroPtmSearchProcess(ZeroPtmMngPtr mng_ptr);

void zeroPtmSearch(SpectrumSetPtr spec_set_ptr, int type,
                   ProteoformPtrVec &form_ptr_vec, ZeroPtmMngPtr mng_ptr,
                   PrSMPtrVec &prsms);

} /* namespace_prot */

#endif
