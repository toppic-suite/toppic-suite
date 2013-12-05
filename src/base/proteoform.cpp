#include <sstream>

#include "base/proteoform.hpp"

#include <log4cxx/logger.h>

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Proteofom"));

Proteoform::Proteoform(std::string name, ResSeqPtr res_seq_ptr) {
  name_ = name;
  res_seq_ptr_ = res_seq_ptr;
  start_pos_ = 0;
  end_pos_ = res_seq_ptr->getLen() - 1;
  LOG4CXX_TRACE(logger, "start bp spec generation");
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

std::string Proteoform::toString() {
  std::stringstream s;
  s<< "Name: " << name_ << std::endl;
  s<< "Begin pos: " << start_pos_ << std::endl;
  s<< "End pos: " << end_pos_ << std::endl;
  s<< "String: " << res_seq_ptr_->toString();
  return s.str();

}

} /* namespace prot */

