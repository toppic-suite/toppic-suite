#include "spec/prm_base_type.hpp"

namespace prot {

const PrmBaseTypePtr PrmBaseType::ORIGINAL = PrmBaseTypePtr(new PrmBaseType("Original"));
const PrmBaseTypePtr PrmBaseType::REVERSED = PrmBaseTypePtr(new PrmBaseType("Reversed"));

}
