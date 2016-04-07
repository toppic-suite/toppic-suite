
#ifndef PROT_TAG_FILTER
#define PROT_TAG_FILTER

#include "base/proteoform.hpp"
#include "spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "tag_filter_mng.hpp"

namespace prot {

class TagFilter {
 public:
  TagFilter(const ProteoformPtrVec &proteo_ptrs,
            TagFilterMngPtr mng_ptr);

  SimplePrsmPtrVec getBestMatch(const PrmMsPtrVec &ms_ptr_vec);

 private:
  TagFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  //CompShiftPtr index_ptr_;

  SimplePrsmPtrVec compute(const PrmMsPtrVec &ms_ptr_vec);
};

typedef std::shared_ptr<TagFilter> TagFilterPtr;
} /* namespace prot */

#endif
