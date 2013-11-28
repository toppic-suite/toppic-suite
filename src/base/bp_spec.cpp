/*
 * bp_spec.cpp
 *
 *  Created on: Nov 26, 2013
 *      Author: xunlikun
 */


#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "log4cxx/logger.h"

#include "bp_spec.hpp"
#include "ion_type.hpp"
#include "mass_comparator.hpp"

namespace prot {

log4cxx::LoggerPtr bs_logger(log4cxx::Logger::getLogger("FastaReader"));

BpSpec::BpSpec(RSPtr res_seq){
	res_seq_=res_seq;
	initBreakPoints();
}

void BpSpec::initBreakPoints(){
	int ext_len= res_seq_->getLen()+1;
	if(ext_len <= 1){
		ext_len = 2;
	}
	break_point_ptr_vec_[0] = BreakPointPtr(new BreakPoint(0,res_seq_->getResMassSum()));
	double prm = 0;
	for(int i=0;i<res_seq_->getLen()-1;i++){
		prm += res_seq_->getResiduePtr(i)->getMass();
		double srm = res_seq_->getResMassSum()-prm;
		if(srm <0){
			LOG4CXX_WARN(bs_logger, "prms is larger than totle mass! ");
		}
		break_point_ptr_vec_.push_back(BreakPointPtr(new BreakPoint(prm,srm)));
	}
	break_point_ptr_vec_.push_back(BreakPointPtr(new BreakPoint(res_seq_->getResMassSum(),0)));
}

std::vector<double> BpSpec::getBreakPointMasses(IonTypePtr ion_type){
	std::vector<double> bpmass_vec;
	if (ion_type->isNTerm()) {
		for (unsigned int i = 0; i < break_point_ptr_vec_.size(); i++) {
			bpmass_vec.push_back(break_point_ptr_vec_[i]->getNTermMass(ion_type));
		}
	}
  else {
		for (unsigned int i = 0; i < break_point_ptr_vec_.size(); i++) {
			bpmass_vec.push_back(break_point_ptr_vec_[i]->getCTermMass(ion_type));
		}
	}
	std::sort(bpmass_vec.begin(),bpmass_vec.end(),MassComparator::up);
	return bpmass_vec;
}

void BpSpec::addBreakPointMass(double mass,double seq_mass,double min_mass,std::vector<double> mass_vec){
	if (mass >= min_mass &&  mass <= seq_mass - min_mass){
		mass_vec.push_back(mass);
	}
}

std::vector<double> BpSpec::getBreakPointMasses(double n_term_shift,double c_term_shift,double min_mass,
		IonTypePtr ion_type_ptr_n,IonTypePtr ion_type_ptr_c){
	std::vector<double> result;
	result.push_back(0.0);
	double new_seq_mass = getResSeq()->getSeqMass()+n_term_shift+c_term_shift;
	//n
	for(unsigned int i=0; i < break_point_ptr_vec_.size();i++){
		double mass = break_point_ptr_vec_[i]->getNTermMass(ion_type_ptr_n)+n_term_shift;
		addBreakPointMass(mass,new_seq_mass,min_mass,result);
	}
	//c
	for(unsigned int i=0; i < break_point_ptr_vec_.size();i++){
		double mass = break_point_ptr_vec_[i]->getCTermMass(ion_type_ptr_c)+n_term_shift;
		addBreakPointMass(mass,new_seq_mass,min_mass,result);
	}
	result.push_back(new_seq_mass);
	std::sort(result.begin(),result.end(),MassComparator::up);

	return result;
}

std::vector<int> BpSpec::getScaledMass(double scale,IonTypePtr ion_type){
	std::vector<int> result;
	if (ion_type->getName().compare("B")==0 || ion_type->getName().compare("C")==0){
		for(unsigned int i=0; i < break_point_ptr_vec_.size();i++){
			double value = break_point_ptr_vec_[i]->getNTermMass(ion_type)*scale;
			result.push_back(std::floor(value+0.5));
		}
	}
	if (ion_type->getName().compare("Y")==0 || ion_type->getName().compare("Z_DOT")==0){
		for(unsigned int i=0; i < break_point_ptr_vec_.size();i++){
			double value = break_point_ptr_vec_[i]->getCTermMass(ion_type)*scale;
			result.push_back(std::floor(value+0.5));
		}
	}
	return result;
}

int getFirstResPos(double n_term_shift,std::vector<double> extbmasses){
	double trunc_len = - n_term_shift;
	int best_pos = -1;
	double best_shift = 2139095040.0;
	for(unsigned int i = 0; i < extbmasses.size();i++){
		if(abs(extbmasses[i] - trunc_len) < best_shift){
			best_pos = i;
			best_shift = abs(extbmasses[i] - trunc_len);
		}
	}
	return best_pos;
}

int getLastResPos(double c_term_shift,std::vector<double> extbmasses){
	double trunc_len = -c_term_shift;
	int best_pos = -1;
	double best_shift = 2139095040.0;
	double pep_mass = extbmasses[extbmasses.size()-1];
	for(unsigned int i=0;i<extbmasses.size();i++){
		if (abs(pep_mass-extbmasses[i]-trunc_len)<best_shift){
			best_pos=i;
			best_shift = abs(pep_mass-extbmasses[i]-trunc_len);
		}
	}
	if(best_pos < 0){
		LOG4CXX_ERROR(bs_logger, "get last residue position error! ");
		throw "get last residue position error!";
	}
	return best_pos - 1;
}

BpSpecPtrVec readBpspecDb(RSPtrVec rs_list){
	BpSpecPtrVec bpspec_ptr_list;
	for(unsigned int i =0; i<rs_list.size();i++){
		bpspec_ptr_list.push_back(BpSpecPtr(new BpSpec(rs_list[i])));
	}
	return bpspec_ptr_list;
}

} /* namespace prot */
