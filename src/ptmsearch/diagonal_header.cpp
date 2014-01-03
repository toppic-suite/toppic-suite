/*
 * diagonal_header.cpp
 *
 *  Created on: Dec 30, 2013
 *      Author: xunlikun
 */

#include <ptmsearch/diagonal_header.hpp>

namespace prot {

DiagonalHeader::DiagonalHeader(double n_term_shift,bool n_strict,bool c_strict,bool n_trunc,bool c_trunc){
	prot_N_term_shift_ =n_term_shift;
	n_strict_=n_strict;
	c_strict_=c_strict;
	n_trunc_ =n_trunc;
	c_trunc_ =c_trunc;
}
DiagonalHeaderPtr DiagonalHeader::clone(){
	DiagonalHeaderPtr cloned = DiagonalHeaderPtr(new DiagonalHeader(prot_N_term_shift_,n_strict_,c_strict_,n_trunc_,c_trunc_));
	cloned->setId(id_);
	cloned->setTruncFirstResPos(trunc_first_res_pos_);
	cloned->setMatchFirstResPos(match_first_res_pos_);
	cloned->setPepNTermShift(pep_N_term_shift_);
	cloned->setProtNTermAllowTrunc(prot_N_term_allow_trunc_);
	cloned->setProtNTermAllowMod(prot_N_term_allow_mod_);
	cloned->setPepNTermAllowMod(pep_N_term_allow_mod_);
	cloned->setAlignPrefix(is_align_prefix_);
	cloned->setAlignsuffix(is_align_suffix_);
	cloned->setTruncLastResPos(trunc_last_res_pos_);
	cloned->setMatchLastResPos(match_last_res_pos_);
	cloned->setProtCTermShift(prot_C_term_shift_);
	cloned->setPepCTermShift(pep_C_term_shift_);
	cloned->setProtCTermAllowTrunc(prot_C_term_allow_trunc_);
	cloned->setProtCTermAllowMod(prot_C_term_allow_mod_);
	cloned->setPepCTermAllowMod(pep_C_term_allow_mod_);
	return cloned;
}
DiagonalHeaderPtr DiagonalHeader::getShift(DiagonalHeaderPtr shift,int bgn,int end){
	DiagonalHeaderPtr new_shift = shift->clone();
	new_shift->setMatchFirstResPos(bgn);
	new_shift->setMatchLastResPos(end);
	return new_shift;
}

DiagonalHeaderPtrVec getNTermShiftListCommon(std::vector<double> best_shifts){
	DiagonalHeaderPtrVec headers;
	for(int i =0;i<best_shifts.size();i++){
		headers.push_back(DiagonalHeaderPtr(new DiagonalHeader(best_shifts[i],true,false,false,false)));
	}
	return headers;
}

DiagonalHeaderPtrVec getNTermShiftListCompLeft(ProteoformPtr seq,PtmMngPtr mng){
	DiagonalHeaderPtrVec extend_n_term_shifts;
	double shift;
	for(int i=0;i<mng->allow_prot_N_mods_.size();i++){
		if(seq->getResSeqPtr()->allowsMod(mng->allow_prot_N_mods_[i]) && mng->allow_prot_N_mods_[i]->getPepShift()==0){
			shift = mng->allow_prot_N_mods_[i]->getProtShift();
			extend_n_term_shifts.push_back(DiagonalHeaderPtr(new DiagonalHeader(shift,true,false,true,false)));
		}
	}
	return extend_n_term_shifts;
}
DiagonalHeaderPtrVec getNTermShiftListCompRight(ProteoformPtr seq,PrmMsPtr ms_six){
	DiagonalHeaderPtrVec extend_n_term_shifts;
	std::vector<double> ms_masses = prot::getMassList(ms_six);
	std::vector<double> seq_masses = seq->getBpSpecPtr()->getBreakPointMasses(IonTypePtr(new IonType("B",true,0)));
	double shift = ms_masses[ms_masses.size()-1] - seq_masses[seq_masses.size()-1];
	extend_n_term_shifts.push_back(DiagonalHeaderPtr(new DiagonalHeader(shift,true,false,true,false)));
	return extend_n_term_shifts;
}

void setPrefixSuffix(DiagonalHeaderPtr header,double c_shift,ProteoformPtr seq,PtmMngPtr mng){
	header->setProtCTermShift(c_shift);
	std::vector<double> seq_b_masses = seq->getBpSpecPtr()->getBreakPointMasses(IonTypePtr(new IonType("B",true,0)));
	double prot_n_term_shift = header->getProtNTermShift();

}

} /* namespace prot */
