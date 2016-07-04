#ifndef PROT_SPEC_RAW_MS_HPP_
#define PROT_SPEC_RAW_MS_HPP_

#include <memory>
#include <vector>

#include "spec/peak.hpp"
#include "spec/ms.hpp"

namespace prot {

typedef std::shared_ptr<Ms<PeakPtr>> RawMsPtr;
typedef std::vector<RawMsPtr> RawMsPtrVec;

} /* namespace prot */

#endif 
