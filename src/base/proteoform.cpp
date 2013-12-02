#include "base/proteoform.hpp"

namespace prot {

Proteoform::Proteoform(std::string name, ResSeqPtr res_seq_ptr) {
  name_ = name;
  res_seq_ptr_ = res_seq_ptr;
  start_pos_ = 0;
  end_pos_ = res_seq_ptr->getLen() - 1;
  bp_spec_ptr_ = BpSpecPtr(new BpSpec(res_seq_ptr));
}

} /* namespace prot */

