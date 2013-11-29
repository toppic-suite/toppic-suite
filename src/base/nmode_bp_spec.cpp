/*
 * nmode_bp_spec.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: xunlikun
 */

#include <nmode_bp_spec.hpp>
#include "bp_spec.hpp"
#include "GLOBAL.hpp"

namespace prot {
NModeBpSpec::NModeBpSpec(std::string name,ResiduePtrVec residues,BpSpecPtr unmode_bpspec_ptr,ProtModPtr n_mode):BpSpec(ResSeqPtr(new ResidueSeq(name,residues))){
	unmode_bpspec_ptr_ = unmode_bpspec_ptr;
	n_mode_ = n_mode;
}

NModeBpSpecPtr getInstance(BpSpecPtr bp_spec,ProtModPtr n_mod){
	PtmPtr acetylation =prot::getPtmPtrByAbbrName(prot::_G_PtmPtrVec,"ptm_list");
	ReqSeqPtr res_seq = bp_spec->getResSeq();
	//TODO:
	//if(!res_seq.allowMode(nMod)){
	//	return nullptr;
	//}
	std::string name = bp_spec->getResSeq()->getName();
	ResiduePtrVec residues = bp_spec->getResSeq()->getResidues();

	//todo:
	if(n_mod == prot::getProtModPtrByName(prot::_G_ProtModPtrVec,"NONE")){
		return NModeBpSpecPtr(new NModeBpSpec(name,residues,bp_spec,n_mod));
	}
	else if (n_mod == prot::getProtModPtrByName(prot::_G_ProtModPtrVec,"NME")){
		ResiduePtrVec mod_res;
		mod_res.assign(++residues.begin(),residues.end());
		return NModeBpSpecPtr(new NModeBpSpec(name,mod_res,bp_spec,n_mod));
	}
	else if (n_mod == prot::getProtModPtrByName(prot::_G_ProtModPtrVec,"ACETYLATION")){
		ResiduePtrVec mod_res;
		mod_res.assign(residues.begin(),residues.end());
		mod_res[0] = prot::getResiduePtrByAcidPtm(prot::_G_ResiduePtrVec,mod_res[0]->getAcidPtr(),acetylation);
		if (mod_res[0] != nullptr){
			prot::addResidue(prot::_G_ResiduePtrVec,mod_res[0]->getAcidPtr(),acetylation);
		}
		return NModeBpSpecPtr(new NModeBpSpec(name,mod_res,bp_spec,n_mod));
	}
	else if (n_mod == prot::getProtModPtrByName(prot::_G_ProtModPtrVec,"NME_ACETYLATION")){
		ResiduePtrVec mod_res;
		mod_res.assign(++residues.begin(),residues.end());
		mod_res[0] = prot::getResiduePtrByAcidPtm(prot::_G_ResiduePtrVec,mod_res[0]->getAcidPtr(),acetylation);
		if (mod_res[0] != nullptr){
			prot::addResidue(prot::_G_ResiduePtrVec,mod_res[0]->getAcidPtr(),acetylation);
		}
		return NModeBpSpecPtr(new NModeBpSpec(name,mod_res,bp_spec,n_mod));
	}
	return nullptr;
}

NModeBpSpecPtrVec getInstance(BpSpecPtrVec bp_spec_list,ProtModPtrVec n_mode_list){
	NModeBpSpecPtrVec nmode_bpspec_ptr_list;
	for(unsigned int i =0;i<bp_spec_list.size();i++){
		for (unsigned int j=0;j<n_mode_list.size();j++){
			NModeBpSpecPtr nmode_bespec_ptr = prot::getInstance(bp_spec_list[i],n_mode_list[j]);
			if (nmode_bespec_ptr != nullptr){
				nmode_bpspec_ptr_list.push_back(nmode_bespec_ptr);
			}
		}
	}
	return nmode_bpspec_ptr_list;
}

} /* namespace prot */
