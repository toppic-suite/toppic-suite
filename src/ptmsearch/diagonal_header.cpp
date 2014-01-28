/*
 * diagonal_header.cpp
 *
 *  Created on: Dec 30, 2013
 *      Author: xunlikun
 */

#include <ptmsearch/diagonal_header.hpp>
#include "base/prot_mod.hpp"
#include "base/change.hpp"

namespace prot {

DiagonalHeader::DiagonalHeader(double n_term_shift,bool n_strict,bool c_strict,
		bool n_trunc,bool c_trunc){
	prot_N_term_shift_ =n_term_shift;
	n_strict_=n_strict;
	c_strict_=c_strict;
	n_trunc_ =n_trunc;
	c_trunc_ =c_trunc;
}
DiagonalHeaderPtr DiagonalHeader::clone(){
	DiagonalHeaderPtr cloned = DiagonalHeaderPtr(
			new DiagonalHeader(prot_N_term_shift_,n_strict_,c_strict_,n_trunc_,c_trunc_));
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
DiagonalHeaderPtr getShift(DiagonalHeaderPtr shift,int bgn,int end){
	DiagonalHeaderPtr new_shift = shift->clone();
	new_shift->setMatchFirstResPos(bgn);
	new_shift->setMatchLastResPos(end);
	return new_shift;
}

DiagonalHeaderPtrVec getNTermShiftListCommon(std::vector<double> best_shifts){
	DiagonalHeaderPtrVec headers;
	for(unsigned int i =0;i<best_shifts.size();i++){
		headers.push_back(DiagonalHeaderPtr(
				new DiagonalHeader(best_shifts[i],true,false,false,false)));
	}
	return headers;
}

DiagonalHeaderPtrVec getNTermShiftListCompLeft(ProteoformPtr seq,PtmMngPtr mng){
	DiagonalHeaderPtrVec extend_n_term_shifts;
	double shift;
	for(unsigned int i=0;i<mng->allow_prot_N_mods_.size();i++){
		if(mng->allow_prot_N_mods_[i]->allowMod(seq->getResSeqPtr()->getResidues())
				&& mng->allow_prot_N_mods_[i]->getPepShift()==0){
			shift = mng->allow_prot_N_mods_[i]->getProtShift();
			extend_n_term_shifts.push_back(DiagonalHeaderPtr(
					new DiagonalHeader(shift,true,false,true,false)));
		}
	}
	return extend_n_term_shifts;
}

DiagonalHeaderPtrVec getNTermShiftListCompRight(ProteoformPtr seq,PrmMsPtr ms_six){
	DiagonalHeaderPtrVec extend_n_term_shifts;
	std::vector<double> ms_masses = prot::getMassList(ms_six);
	std::vector<double> seq_masses = seq->getBpSpecPtr()->getBreakPointMasses(
			IonTypePtr(new IonType("B",true,0)));
	double shift = ms_masses[ms_masses.size()-1] - seq_masses[seq_masses.size()-1];
	extend_n_term_shifts.push_back(DiagonalHeaderPtr(
			new DiagonalHeader(shift,false,true,false,true)));
	return extend_n_term_shifts;
}

void setPrefixSuffix(DiagonalHeaderPtr &header,double c_shift,ProteoformPtr seq,
		PtmMngPtr mng){
	header->setProtCTermShift(c_shift);

	std::vector<double> seq_b_masses = seq->getBpSpecPtr()->getBreakPointMasses(
			IonTypePtr(new IonType("B",true,0)));
	double prot_n_term_shift = header->getProtNTermShift();
	double prot_c_term_shift = header->getProtCTermShift();
	int trunc_first_res_pos = prot::getFirstResPos(prot_n_term_shift,seq_b_masses);

	header->setTruncFirstResPos(trunc_first_res_pos);
	int trunc_last_res_pos = prot::getLastResPos(prot_c_term_shift,seq_b_masses);
//	std::cout<<trunc_last_res_pos<<std::endl;
	header->setTruncLastResPos(trunc_last_res_pos);
	double pep_n_term_shift = prot_n_term_shift + seq_b_masses[trunc_first_res_pos];
	header->setPepNTermShift(pep_n_term_shift);
	double pep_c_term_shift = prot_c_term_shift + seq_b_masses[seq_b_masses.size()-1]
	                                 -seq_b_masses[trunc_last_res_pos+1];
	header->setPepCTermShift(pep_c_term_shift);

	prot::setProtTermMod(header,seq,mng);
	prot::setProtTermTrunc(header,seq,mng);
	prot::setAlignPrefSuffic(header,mng);
	prot::setPepTermMode(header,mng);
}

void setProtTermMod(DiagonalHeaderPtr &header,ProteoformPtr seq,PtmMngPtr mng){
	int trunc_len = header->getTruncFirstResPos();

	ResSeqPtr resseq= seq->getResSeqPtr();

	ProtModPtr mod;
	if(header->isNTrunc()){
		mod = findProtTermMod(mng->allow_prot_N_mods_,
				trunc_len,resseq,
				header->getPepNTermShift(),
				mng->test_term_mod_error_toerance_);
	}
	header->setProtNTermAllowMod(mod);
	trunc_len = seq->getResSeqPtr()->getLen()-1 -header->getTruncLastResPos();
	resseq = resseq->getSubResidueSeq(header->getTruncLastResPos()+1,resseq->getLen()-1);
	mod = nullptr;
	if(header->isCTrunc()){
		mod = findProtTermMod(mng->allow_prot_C_mods_,
				trunc_len,resseq,
				header->getPepCTermShift(),
				mng->test_term_mod_error_toerance_);
	}
	header->setProtCTermAllowMod(mod);
}

void setProtTermTrunc(DiagonalHeaderPtr &header,ProteoformPtr seq,PtmMngPtr mng){
	TruncPtr trunc = prot::findProtNTermTrunc(
			seq->getResSeqPtr(),
			header->getTruncFirstResPos(),
			mng->allow_prot_N_truncs_);
	header->setProtNTermAllowTrunc(trunc);
	trunc = prot::findProtCTermTrunc(
			seq->getResSeqPtr(),
			header->getTruncLastResPos(),
			mng->allow_prot_C_truncs_);
	header->setProtCTermAllowTrunc(trunc);
}

void setPepTermMode(DiagonalHeaderPtr &header,PtmMngPtr mng){
	PtmPtr mod;
	if(header->isNTrunc()){
		mod = findPepTermMod(
				mng->allow_pep_N_mods_,
				header->getPepNTermShift(),
				mng->test_term_mod_error_toerance_);
	}
	header->setPepNTermAllowMod(mod);
	mod =nullptr;
	if(header->isCTrunc()){
		mod = findPepTermMod(mng->allow_pep_C_mods_,
				header->getPepCTermShift(),
				mng->test_term_mod_error_toerance_);
	}
	header->setPepCTermAllowMod(mod);
}
ProtModPtr findProtTermMod(ProtModPtrVec mods,int trunc_len,
		ResSeqPtr res_seq,double pep_term_shift,double tolerance){
	for(unsigned int i=0;i<mods.size();i++){
		if(mods[i]->getTruncPtr()->isSameTrunc(trunc_len,res_seq)
				&& std::abs(mods[i]->getPepShift()-pep_term_shift<=tolerance)){
			return mods[i];
		}
	}
	return nullptr;
}

PtmPtr findPepTermMod(PtmPtrVec mods,double shift,double tolerance){
	for(unsigned int i=0;i<mods.size();i++){
		if(shift>=mods[i]->getMonoMass() && shift-mods[i]->getMonoMass()<= tolerance){
			return mods[i];
		}
		if(shift<mods[i]->getMonoMass() && mods[i]->getMonoMass()-shift<= tolerance){
			return mods[i];
		}
	}
	return nullptr;
}

void setAlignPrefSuffic(DiagonalHeaderPtr &header,PtmMngPtr mng){
	bool is_prefix = prot::isAlignPrefix(header->getProtNTermAllowTrunc(),
			header->getPepNTermShift(),mng->prefix_suffix_shift_thesh_);
	header->setAlignPrefix(is_prefix);
	bool is_suffix = prot::isAlignSuffix(header->getProtCTermAllowTrunc(),
			header->getPepCTermShift(),mng->prefix_suffix_shift_thesh_);
	header ->setAlignsuffix(is_suffix);
}

DiagonalHeaderPtrVec getNTermShiftListTruncPrefix(ProteoformPtr seq){
	DiagonalHeaderPtrVec extend_n_term_shift;
	std::vector<double> seq_masses = seq->getBpSpecPtr()->getBreakPointMasses(
			IonTypePtr(new IonType("B",true,0)));
	double shift;

	for(unsigned int i=1;i<seq_masses.size();i++){
		shift = - seq_masses[i];
		extend_n_term_shift.push_back(
				DiagonalHeaderPtr(new DiagonalHeader(shift,true,false,true,false)));
	}

	return extend_n_term_shift;
}
DiagonalHeaderPtrVec getNTermShiftListTruncsuffix(PrmMsPtr ms,ProteoformPtr seq){
	DiagonalHeaderPtrVec extend_n_term_shift;
	std::vector<double> ms_masses = prot::getMassList(ms);
	std::vector<double> seq_masses = seq->getBpSpecPtr()->getBreakPointMasses(
			IonTypePtr(new IonType("B",true,0)));
	double shift;

	for(unsigned int i=1;i<seq_masses.size();i++){
		shift = ms_masses[ms_masses.size() - 1]-seq_masses[i];
		extend_n_term_shift.push_back(
				DiagonalHeaderPtr(new DiagonalHeader(shift,true,false,true,false)));
	}
	return extend_n_term_shift;
}
DiagonalHeaderPtrVec get1dHeaders(DiagonalHeaderPtrVec2D headers){
	DiagonalHeaderPtrVec header_list;
		for(unsigned int i =0;i<headers.size();i++){
			for(unsigned int j=0;j<headers[i].size();j++){
				header_list.push_back(headers[i][j]);
			}
		}
	return header_list;
}

bool getNAcetylation(DiagonalHeaderPtrVec headers){
	if(headers.size()==0){
		return false;
	}
	ProtModPtr mod = headers[0]->getProtNTermAllowMod();
//	std::cout<<mod<<std::endl;
	if(mod==nullptr||mod->getName().compare("ACETYLATION")==0
			||mod->getName().compare("NME_ACETYLATION")==0){
		return true;
	}
	return false;
}

ChangePtrVec getChanges(DiagonalHeaderPtrVec headers,int first,int last,PtmPtrVec ptm_list){

	ChangePtrVec change_list;
	if(headers[0]->getPepNTermAllowMod() != PtmFactory::findEmptyPtmPtr()){
		if(getNAcetylation(headers)){
			ProtModPtr aptm = headers[0]->getProtNTermAllowMod();
			change_list.push_back(ChangePtr(
					new Change(first,headers[0]->getMatchFirstResPos(),
							PROTEIN_VARIABLE_CHANGE,
							headers[0]->getPepNTermShift(),
							aptm==nullptr?nullptr:aptm->getPtmPtr())));
		}
		else{
			change_list.push_back(ChangePtr(
					new Change(first,headers[0]->getMatchFirstResPos(),
							UNEXPECTED_CHANGE,
							headers[0]->getPepNTermShift(),
							nullptr)));
		}
	}
	for(unsigned int i =0;i<headers.size()-1;i++){
		change_list.push_back(ChangePtr(
				new Change(first,headers[0]->getMatchFirstResPos(),
						UNEXPECTED_CHANGE,
						headers[0]->getPepNTermShift(),
						nullptr)));
	}
	DiagonalHeaderPtr lastHeader = headers[headers.size()-1];
	if(lastHeader->getPepCTermAllowMod() != PtmFactory::findEmptyPtmPtr()){
		change_list.push_back(ChangePtr(
				new Change(lastHeader->getMatchLastResPos()+1,
						last+1,UNEXPECTED_CHANGE,
						lastHeader->getPepCTermShift(),
						nullptr)));
	}

	return change_list;
}

} /* namespace prot */
