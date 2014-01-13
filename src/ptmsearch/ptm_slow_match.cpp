/*
 * ptm_slow_match.cpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#include <ptmsearch/ptm_slow_match.hpp>

namespace prot {

PtmSlowMatch::PtmSlowMatch(ProteoformPtr seq,SpectrumSetPtr spectrum_set,CompShiftLowMemPtr comp_shift,PtmMngPtr mng){
	mng_=mng;
	deconv_ms_ = spectrum_set->getDeconvMs();
	ms_six_=spectrum_set->getSpSix();
	ms_three_ = spectrum_set->getSpThree();
	seq_=seq;
	comp(comp_shift);
}

void PtmSlowMatch::comp(CompShiftLowMemPtr comp_shift){
	double scale = mng_->ptm_fast_filter_scale_;
	std::vector<std::vector<int>> sp_masses = prot::getIntMassErrorList(ms_six_,scale,true,false);
	std::vector<double> best_shift= comp_shift->findBestShift(sp_masses[0],sp_masses[1],seq_->getBpSpecPtr()->getScaledMass(scale,IonTypePtr(new IonType("B",true,0))),mng_->n_top_diagonals_,mng_->min_diagonal_gap_,scale);
	DiagonalHeaderPtrVec n_term_shifts =getNTermShiftList(best_shift,ms_six_,seq_,mng_);
	BasicDiagPairDiagPtrVec diagonals = prot::getDiagonals(n_term_shifts,ms_six_,seq_,mng_);
	std::vector<double> ms_masses = prot::getMassList(ms_six_);
	std::vector<double> seq_masses = seq_->getBpSpecPtr()->getBreakPointMasses(IonTypePtr(new IonType("B",true,0)));
	PSAlignPtr align = PSAlignPtr(new PSAlign(ms_masses,seq_masses,diagonals,mng_));
	std::vector<int> types = {SEMI_ALIGN_TYPE_COMPLETE,SEMI_ALIGN_TYPE_PREFIX,SEMI_ALIGN_TYPE_SUFFIX,SEMI_ALIGN_TYPE_INTERNAL};

//	std::vector<std::vector<double>> result_scores;
//	DiagonalHeaderPtrVec3D result_headers;

	for(int i=0;i<types.size();i++){
		align->compute(types[i]);
		result_scores_.push_back(align->getAlignScr());
		result_headers_.push_back(align->getResult());
	}

	for(int i=0;i<4;i++){
		std::vector<double> temp;
		for(int j=0;j<=mng_->n_unknown_shift_;j++){
			temp.push_back(0);
		}
		result_deltas_.push_back(temp);
	}
}
double PtmSlowMatch::getScr(int shiftnum,int type){
	return result_scores_[type][shiftnum];
}

PrSMPtr PtmSlowMatch::geneResult(int shift_num,int type){
	DiagonalHeaderPtrVec headers=result_headers_[type][shift_num];
	double refine_prec_mass = ms_three_->getHeaderPtr()->getPrecMonoMass()+result_deltas_[type][shift_num];
	int first_pos = headers[0]->getTruncFirstResPos();
	int last_pos = headers[headers.size()-1]->getTruncLastResPos();
	DiagonalHeaderPtrVec refined_headers = prot::refineHeadersBgnEnd(first_pos,seq_,deconv_ms_,ms_three_,mng_,headers);
	if(refined_headers.size()==0){
		return nullptr;
	}

	ChangePtrVec changes = prot::getChanges(headers,first_pos,last_pos);
//	Proteoform(DbResSeqPtr db_res_seq_ptr, ProtModPtr prot_mod_ptr,
//	             ResSeqPtr res_seq_ptr, int start_pos, int end_pos,
//	             ChangePtrVec change_list);


	ProteoformPtr protein = ProteoformPtr(new Proteoform(seq_->getDbResSeqPtr(),seq_->getProtModPtr(),seq_->getResSeqPtr(),first_pos,last_pos,changes) );
	return PrSMPtr(new PrSM(protein,deconv_ms_,refine_prec_mass,0,mng_->sp_para_));
}

DiagonalHeaderPtrVec PtmSlowMatch::getNTermShiftList(std::vector<double> best_shift,PrmMsPtr ms_six,ProteoformPtr seq,PtmMngPtr mng){
	DiagonalHeaderPtrVec headers = prot::getNTermShiftListCommon(best_shift);
	DiagonalHeaderPtrVec n_term_shifts_comp_left = prot::getNTermShiftListCompLeft(seq,mng);
	DiagonalHeaderPtrVec n_term_shifts_comp_right = prot::getNTermShiftListCompRight(seq,ms_six);
	DiagonalHeaderPtrVec extend_n_term_shifts;
	std::vector<double> ms_masses = prot::getMassList(ms_six);
	std::vector<double> seq_masses = seq->getBpSpecPtr()->getBreakPointMasses(IonTypePtr(new IonType("B",true,0)));
	double shift;

	for(int i=1;i<seq_masses.size();i++){
		shift = - seq_masses[i];
		if(found(shift,headers,mng)){
			extend_n_term_shifts.push_back(DiagonalHeaderPtr(new DiagonalHeader(shift,true,false,true,false)));
		}
	}

	for(int i=1;i<seq_masses.size();i++){
		shift = ms_masses[ms_masses.size()-1] - seq_masses[i];
		if(found(shift,headers,mng)){
			extend_n_term_shifts.push_back(DiagonalHeaderPtr(new DiagonalHeader(shift,true,false,true,false)));
		}
	}

	for(int i=0;i<extend_n_term_shifts.size();i++){
		headers.push_back(extend_n_term_shifts[i]);
	}

	return headers;
}

bool PtmSlowMatch::found(double shift,DiagonalHeaderPtrVec headerlist,PtmMngPtr mng){
	for(int i=0;i<headerlist.size();i++){
		if(std::abs(shift-headerlist[i]->getProtNTermShift())<= mng->extend_diagonal_error_tolerance_){
			return true;
		}
	}
	return false;
}
} /* namespace prot */
