#ifndef PROT_TRUNC_UTIL_H_
#define PROT_TRUNC_UTIL_H_

#include <vector>

#include "trunc.hpp"

namespace prot {


TruncPtrVec getTruncPtrVecInstance(AcidPtrVec &acid_ptr_vec, 
                                         const char* file_name);

TruncPtr getTruncPtrByName(TruncPtrVec &trunc_ptr_vec, 
                         const std::string &name);
}
#endif
