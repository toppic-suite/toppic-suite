#include "spec/prm_break_type.hpp"

namespace prot {

const PrmBreakTypePtr PrmBreakType::NONE = PrmBreakTypePtr(new PrmBreakType("NONE"));
const PrmBreakTypePtr PrmBreakType::N_TERM = PrmBreakTypePtr(new PrmBreakType("N_term"));
const PrmBreakTypePtr PrmBreakType::C_TERM = PrmBreakTypePtr(new PrmBreakType("C_term"));
const PrmBreakTypePtr PrmBreakType::BOTH = PrmBreakTypePtr(new PrmBreakType("Both"));

}
