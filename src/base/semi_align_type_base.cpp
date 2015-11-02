#include "base/semi_align_type_base.hpp"

namespace prot {

SemiAlignTypePtr SemiAlignTypeBase::semi_align_type_complete_ 
    = SemiAlignTypePtr(new SemiAlignType("COMPLETE", 0));
SemiAlignTypePtr SemiAlignTypeBase::semi_align_type_prefix_
    = SemiAlignTypePtr(new SemiAlignType("PREFIX", 1));
SemiAlignTypePtr SemiAlignTypeBase::semi_align_type_suffix_  
    = SemiAlignTypePtr(new SemiAlignType("SUFFIX", 2));
SemiAlignTypePtr SemiAlignTypeBase::semi_align_type_internal_
    = SemiAlignTypePtr(new SemiAlignType("INTERNAL", 3));

}
