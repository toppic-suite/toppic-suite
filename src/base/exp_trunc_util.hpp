#ifndef PROT_EXP_TRUNC_UTIL_H_
#define PROT_EXP_TRUNC_UTIL_H_

#include <vector>

#include "exp_trunc.hpp"

namespace prot {


ExpTruncPtrVec getExpTruncPtrVecInstance(AcidPtrVec &acid_ptr_vec, 
                                         const char* file_name);

ExpTruncPtr getExpTruncPtrByName(ExpTruncPtrVec &exp_trunc_ptr_vec, 
                         const std::string &name);
}
#endif
