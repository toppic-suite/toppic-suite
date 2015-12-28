#ifndef PROT_ONE_PTM_SEARCH_HPP_
#define PROT_ONE_PTM_SEARCH_HPP_

#include <array>

#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "oneptmsearch/one_ptm_search_mng.hpp"

namespace prot {

class OnePtmSearch {
 public:
  static void process(OnePtmSearchMngPtr mng_ptr);
};

} /* namespace_prot */

#endif
