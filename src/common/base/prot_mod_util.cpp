//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <string>

#include "common/util/logger.hpp"
#include "common/base/mod_base.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/trunc_util.hpp"
#include "common/base/prot_mod_util.hpp"

namespace toppic {

namespace prot_mod_util {

bool allowMod(ProtModPtr prot_mod_ptr, const ResiduePtrVec &residues) {
  // Case 1. no protein modification
  if (prot_mod_ptr == ProtModBase::getProtModPtr_NONE()) {
    return true;
  } 
  
  // Case 2. N-terminal methionine acetylation
  if (prot_mod_ptr == ProtModBase::getProtModPtr_M_ACETYLATION()) {
    int mod_pos = prot_mod_ptr->getModPos();
    if (mod_pos >= static_cast<int>(residues.size())) {
      // LOG_DEBUG("pos false");
      return false;
    }
    ModPtr mod_ptr = prot_mod_ptr->getModPtr();
    if (residues[mod_pos] != mod_ptr->getOriResiduePtr()) {
      // LOG_DEBUG("mod false");
      return false;
    }
    return true;
  } 

  // Case 3. NME and NME acetylation
  // check truncation
  if (!trunc_util::isValidTrunc(prot_mod_ptr->getTruncPtr(), residues)) {
    return false;
  }
  ModPtr mod_ptr = prot_mod_ptr->getModPtr();
  if (mod_ptr != ModBase::getNoneModPtr()) {
    // if NME_acetylation
    int mod_pos = prot_mod_ptr->getModPos();
    if (mod_pos >= static_cast<int>(residues.size())) {
      return false;
    }
    if (residues[mod_pos] != mod_ptr->getOriResiduePtr()) {
      return false;
    }
  }
  return true;
}

} // namespace toppic_mod_util

}  // namespace toppic
