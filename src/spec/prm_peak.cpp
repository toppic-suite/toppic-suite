/*
 * prm_peak.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#include <spec/prm_peak.hpp>

namespace prot {
PrmPeak::PrmPeak(DeconvPeakPtr base_peak,PrmPeakTypePtr base_type,double mono_mass,double score):Peak(mono_mass,base_peak->getIntensity()){
	base_peak_=base_peak;
	base_type_=base_type;
	mono_mass_=mono_mass;
	score_=score;
	strict_tolerance_=0;
	n_strict_c_relax_tolerance_=0;
	n_relax_c_strict_tolerance_=0;
}

void PrmPeak::addNghbEdge(DeconvPeakPtr peak,double offset,SupportPeakTypePtr peak_type,double score){
	score_ +=score;
	neighbor_list_.push_back(SupportPeakPtr(new SupportPeak(peak,offset,score,peak_type)));
}

int PrmPeak::getBreakType(SupportPeakTypePtrVec support_peak_type_list){
	int break_type = 0;
	for(unsigned int i=0;i<neighbor_list_.size();i++){
		if(neighbor_list_[i]->getPeakType() == prot::getSupportPeakTypePtrByName(support_peak_type_list,"N_TERM")){
			if(break_type == 0){
				break_type = 1;
			}
			else if(break_type == 2){
				break_type = 3;
			}
		}
		else{
			if(break_type == 0){
				break_type = 2;
			}
			else if(break_type == 1){
				break_type = 3;
			}
		}
	}
	return break_type;
}
} /* namespace prot */
