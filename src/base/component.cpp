/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */

#include "component.hpp"

namespace prot {

Component::Component (const char* acid_file_name,
                      const char* ptm_file_name,
                      const char* residue_file_name) {
  acid_list_ = getAcidPtrVecInstance(acid_file_name);
  ptm_list_ = getPtmPtrVecInstance(ptm_file_name);
  residue_list_ = getResiduePtrVecInstance(acid_list_,
                                         ptm_list_,
                                         residue_file_name);
}

}

