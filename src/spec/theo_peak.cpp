/*
 * theo_peak.cpp
 *
 *  Created on: Nov 28, 2013
 *      Author: xunlikun
 */

#include <algorithm>
#include "base/proteoform.hpp"
#include "spec/theo_peak.hpp"

namespace prot {

TheoPeak::TheoPeak(IonPtr ion,double unmod_mass,double shift):
  Peak(unmod_mass + shift, 1.0) {
	ion_ = ion;
	unmod_mass_ = unmod_mass;
	shift_ = shift;
}

TheoPeakPtrVec getTheoPeak(BpSpecPtr pep,ActivationPtr type, NeutralLossPtr neutral_loss_ptr, 
                           double n_term_shift,double c_term_shift,int bgn,int end,double min_mass){
	TheoPeakPtrVec theo_peaks;
	BreakPointPtrVec bps = pep->getBreakPointPtrVec();
	double new_seq_mass = pep->getRSMass()+n_term_shift;
	IonTypePtr n_ion_type = type->getNIonType();
	for(int i =bgn;i<=end;i++){
		double n_mass = bps[i]->getNTermMass(n_ion_type);
		double new_mass = n_mass + n_term_shift;
		if(new_mass >= min_mass && new_mass <= new_seq_mass - min_mass){
			IonPtr ion = IonPtr(new Ion(0,i,i,n_ion_type, neutral_loss_ptr));
			theo_peaks.push_back(TheoPeakPtr(new TheoPeak(ion,n_mass,n_term_shift)));
		}
	}

	IonTypePtr c_ion_type = type->getCIonType();
	for(int i =bgn;i<=end;i++){
		double c_mass = bps[i]->getCTermMass(c_ion_type);
		double new_mass = c_mass + c_term_shift;
		if(new_mass >= min_mass && new_mass <= new_seq_mass - min_mass){
			IonPtr ion = IonPtr(new Ion(0,i,bps.size()-i-1,c_ion_type,neutral_loss_ptr));
			theo_peaks.push_back(TheoPeakPtr(new TheoPeak(ion,c_mass,c_term_shift)));
		}
	}
	std::sort(theo_peaks.begin(),theo_peaks.end(),theo_peak_down);
	return theo_peaks;
}

TheoPeakPtrVec getProteoformTheoPeak(ProteoformPtr proteoform_ptr, 
                                     ActivationPtr activation_ptr,
                                     double min_mass) {

  BpSpecPtr bp_ptr = proteoform_ptr->getBpSpecPtr();
  TheoPeakPtrVec all_peaks;
  SegmentPtrVec segments = proteoform_ptr->getSegmentPtrVec();
  for (unsigned int i = 0; i < segments.size(); i++) {
    TheoPeakPtrVec  peaks = getTheoPeak(bp_ptr, activation_ptr, NeutralLossPtr(nullptr), 
                                        segments[i]->getPepNTermShift(),
                                        segments[i]->getPepCTermShift(), segments[i]->getLeftBpPos(),
                                        segments[i]->getRightBpPos(), min_mass);
     all_peaks.insert(all_peaks.end(), peaks.begin(), peaks.end());
  }
  return all_peaks;
}

} /* namespace prot */
