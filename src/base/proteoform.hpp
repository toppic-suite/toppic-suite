#ifndef PROT_PROTEOFORM_HPP_
#define PROT_PROTEOFORM_HPP_

#include "residue_seq.hpp"
#include "bp_spec.hpp"
#include "change.hpp"

namespace prot {

class Proteoform {
public:
	Proteoform(ResSeqPtr res_seq_ptr);

	ResSeqPtr getResSeqPtr() {return res_seq_ptr_;}

	BpSpecPtr getBpSpecPtr() {return bp_spec_ptr_;}

  int getStartPos() {return start_pos_;}

  int getEndPos() {return end_pos_;}

  std::vector<Change> getChangeList() {return change_list_;}

private:
	ResSeqPtr res_seq_ptr_;

  BpSpecPtr bp_spec_ptr_;

  int start_pos_;

  int end_pos_;

  std::vector<Change> change_list_;
};

} /* namespace prot */

#endif /* BP_SPEC_HPP_ */
