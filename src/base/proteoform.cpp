#include "base/proteoform.hpp"

namespace prot {

Proteoform::Proteoform(std::string name, ResSeqPtr res_seq_ptr) {
  name_ = name;
  res_seq_ptr_ = res_seq_ptr;
  start_pos_ = 0;
  end_pos_ = res_seq_ptr->getLen() - 1;
  bp_spec_ptr_ = BpSpecPtr(new BpSpec(res_seq_ptr));
  for (int i = 0; i < res_seq_ptr->getLen(); i++) {
    PtmPtr ptm_ptr = res_seq_ptr->getResiduePtr(i)->getPtmPtr();
    if (!ptm_ptr->isEmpty()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(i, i+1, FIXED_CHANGE, ptm_ptr->getMonoMass(), ptm_ptr));
      change_list_.push_back(change_ptr);
    }
  }
}

} /* namespace prot */

