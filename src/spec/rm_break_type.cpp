#include "spec/rm_break_type.hpp"

namespace prot {

const RmBreakTypePtr RmBreakType::NONE = RmBreakTypePtr(new RmBreakType("NONE"));
const RmBreakTypePtr RmBreakType::N_TERM = RmBreakTypePtr(new RmBreakType("N_term"));
const RmBreakTypePtr RmBreakType::C_TERM = RmBreakTypePtr(new RmBreakType("C_term"));
const RmBreakTypePtr RmBreakType::BOTH = RmBreakTypePtr(new RmBreakType("Both"));

}
