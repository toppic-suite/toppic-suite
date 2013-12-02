#ifndef PROT_PROTEOFORM_HPP_
#define PROT_PROTEOFORM_HPP_

#include "base/residue_seq.hpp"
#include "base/bp_spec.hpp"
#include "base/change.hpp"

namespace prot {

class Proteoform {
public:
	Proteoform(std::string name, ResSeqPtr res_seq_ptr);

	ResSeqPtr getResSeqPtr() {return res_seq_ptr_;}

	BpSpecPtr getBpSpecPtr() {return bp_spec_ptr_;}

  int getStartPos() {return start_pos_;}

  int getEndPos() {return end_pos_;}

  std::vector<Change> getChangeList() {return change_list_;}

private:
  std::string name_;

	ResSeqPtr res_seq_ptr_;

  BpSpecPtr bp_spec_ptr_;

  int start_pos_;

  int end_pos_;

  std::vector<Change> change_list_;
};

typedef std::shared_ptr<Proteoform> ProteoformPtr;
typedef std::vector<ProteoformPtr> ProteoformPtrVec;

} /* namespace prot */

#endif /* BP_SPEC_HPP_ */
