#include "spec/base_peak_type.hpp"

namespace prot {

const BasePeakTypePtr BasePeakType::ORIGINAL = BasePeakTypePtr(new BasePeakType("Original"));
const BasePeakTypePtr BasePeakType::REVERSED = BasePeakTypePtr(new BasePeakType("Reversed"));

}
