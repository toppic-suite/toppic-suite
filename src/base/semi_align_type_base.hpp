#ifndef PROT_BASE_SEMI_ALIGN_TYPE_BASE_HPP_
#define PROT_BASE_SEMI_ALIGN_TYPE_BASE_HPP_

#include "base/semi_align_type.hpp"

namespace prot {

class SemiAlignTypeBase {
 private:
  static SemiAlignTypePtr semi_align_type_complete_;
  static SemiAlignTypePtr semi_align_type_prefix_;
  static SemiAlignTypePtr semi_align_type_suffix_;
  static SemiAlignTypePtr semi_align_type_internal_;

 public:
  static SemiAlignTypePtr getCompletePtr() {return semi_align_type_complete_;}
  static SemiAlignTypePtr getPrefixPtr() {return semi_align_type_prefix_;}
  static SemiAlignTypePtr getSuffixPtr() {return semi_align_type_suffix_;}
  static SemiAlignTypePtr getInternalPtr() {return semi_align_type_internal_;}
};

}

#endif
