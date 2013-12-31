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

} /* namespace prot */
