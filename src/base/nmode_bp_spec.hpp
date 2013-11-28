/*
 * nmode_bp_spec.hpp
 *
 *  Created on: Nov 27, 2013
 *      Author: xunlikun
 */

#ifndef NMODE_BP_SPEC_HPP_
#define NMODE_BP_SPEC_HPP_

#include "bp_spec.hpp"
#include "prot_mod.hpp"

namespace prot {

class NModeBpSpec : BpSpec {
public:
	NModeBpSpec(std::string name,ResiduePtrVec residues,BpSpecPtr unmode_bpspec_ptr,ProtModPtr n_mode);
	BpSpecPtr getUnModeBpSpec(){return unmode_bpspec_ptr_;}
	ProtModPtr getNMode(){return n_mode_;}
private:
	BpSpecPtr unmode_bpspec_ptr_;
	ProtModPtr n_mode_;
};

typedef std::shared_ptr<NModeBpSpec> NModeBpSpecPtr;
typedef std::vector<NModeBpSpecPtr> NModeBpSpecPtrVec;

NModeBpSpecPtr getInstance(BpSpecPtr bp_spec,ProtModPtr n_mode);
NModeBpSpecPtrVec getInstance(BpSpecPtrVec bp_spec_list,ProtModPtrVec n_mode_list);

} /* namespace prot */

#endif /* NMODE_BP_SPEC_HPP_ */
