#ifndef PROT_PROT_MOD_UTIL_H_
#define PROT_PROT_MOD_UTIL_H_

#include "prot_mod.hpp"

namespace prot {

ProtModPtrVec getProtModPtrVecInstance(AcidPtrVec &acid_ptr_vec, 
                                       PtmPtrVec &ptm_ptr_vec,
                                       TruncPtrVec &trunc_ptr_vec,
                                       const char* file_name);

ProtModPtr getProtModPtrByName(ProtModPtrVec &prot_mod_ptr_vec, 
                         const std::string &name);
}
#endif
