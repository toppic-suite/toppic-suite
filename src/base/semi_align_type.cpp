#include "base/semi_align_type.hpp"

namespace prot {

SemiAlignType::SemiAlignType(const std::string &name, int id) {
  name_ = name;
  id_ = id;
}

SemiAlignTypePtr SemiAlignTypeFactory::semi_align_type_complete_ 
    = SemiAlignTypePtr(new SemiAlignType("COMPLETE", 0));
SemiAlignTypePtr SemiAlignTypeFactory::semi_align_type_prefix_
    = SemiAlignTypePtr(new SemiAlignType("PREFIX", 1));
SemiAlignTypePtr SemiAlignTypeFactory::semi_align_type_suffix_  
    = SemiAlignTypePtr(new SemiAlignType("SUFFIX", 2));
SemiAlignTypePtr SemiAlignTypeFactory::semi_align_type_internal_
    = SemiAlignTypePtr(new SemiAlignType("INTERNAL", 3));

}
