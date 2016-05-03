#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_util.hpp"

namespace prot {

SimplePrsmPtrVec SimplePrsmUtil::getUniqueMatches(SimplePrsmPtrVec &match_ptrs) {
  std::sort(match_ptrs.begin(), match_ptrs.end(),SimplePrsm::cmpNameIncScoreDec);
  SimplePrsmPtrVec unique_match_ptrs;
  std::string prev_name = "";
  for(size_t i=0;i< match_ptrs.size();i++){
    std::string cur_name = match_ptrs[i]->getSeqName();
    if (cur_name != prev_name) {
      unique_match_ptrs.push_back(match_ptrs[i]);
      prev_name = cur_name;
    }
  }

  return unique_match_ptrs;
}


} /* namespace prot */
