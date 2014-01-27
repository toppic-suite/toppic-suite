/*
 * pair.cpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */
//#include <algorithm>
#include <limits>

#include "prsm/pair.hpp"

namespace prot {

std::vector<double> compPpoDeviation(std::vector<double> ms_masses,
                                     std::vector<double> theo_masses,
                                     double ppo){
	std::vector<double> min_distances;
	for(unsigned int i =0;i< ms_masses.size();i++){
		min_distances.push_back(std::numeric_limits<double>::infinity());
	}
	unsigned int i=0;
	unsigned int j=0;
	while(i<ms_masses.size() && j < theo_masses.size()){
		double d = ms_masses[i] - theo_masses[j];
		if(std::abs(d) < std::abs(min_distances[i])){
			min_distances[i] = d;
		}
		double tolerance = ms_masses[i]*ppo;
		if(increaseIJ(i,j,d,tolerance,ms_masses,theo_masses)){
			i++;
		}
		else{
			j++;
		}
	}

	std::vector<double> min_ppos;

	for(unsigned int i=0;i<ms_masses.size();i++){
		min_ppos.push_back(std::numeric_limits<double>::infinity());
		if(ms_masses[i]>0){
			min_ppos[i] = min_distances[i]/ms_masses[i];
		}
	}
	return min_ppos;
}

double compUniqScore(std::vector<double> ms_masses,
                     std::vector<double> theo_masses,double ppo){
	std::vector<double> min_distances;
	for(unsigned int i =0;i< ms_masses.size();i++){
		min_distances.push_back(std::numeric_limits<double>::infinity());
	}
	unsigned int i=0;
	unsigned int j=0;
	while(i<ms_masses.size()&& j < theo_masses.size()){
		double d = ms_masses[i] - theo_masses[j];
		if(std::abs(d) < std::abs(min_distances[i])){
			min_distances[i] = d;
		}
		double tolerance = ms_masses[i]*ppo;
		if(increaseIJ(i,j,d,tolerance,ms_masses,theo_masses)){
			i++;
		}else{
			j++;
		}
	}

	double score =0;

	for(unsigned int i=0;i<theo_masses.size();i++){
		if(theo_masses[i]>0){
			double error = min_distances[i]/theo_masses[i];
			if(std::abs(error)<=ppo){
				score+=1.0;
			}
		}
	}

	return score;
}

std::vector<double> compPpoDeviation(ExtendMsPtr ms,
                                     TheoPeakPtrVec ions,double ppo){
	std::vector<double> ion_masses;
	for(unsigned int i=0;i<ions.size();i++){
		ion_masses.push_back(ions[i]->getModMass());
	}
	std::vector<double> ms_masses;
	for(unsigned int i=0;i<ms->size();i++){
		ion_masses.push_back(ms->getPeakPtr(i)->getPosition());
	}
	return compPpoDeviation(ms_masses,ion_masses,ppo);
}

double compIonScore(ExtendMsPtr ms,TheoPeakPtrVec ions,
                    double recal,double ppo){
	std::vector<double> ion_masses;
	for(unsigned int i=0;i<ions.size();i++){
		ion_masses.push_back(ions[i]->getModMass());
	}
	std::vector<double> ms_masses;
	for(unsigned int i=0;i<ms->size();i++){
		ms_masses.push_back(ms->getPeakPtr(i)->getPosition()*(1+recal));
	}
	return compUniqScore(ms_masses,ion_masses,ppo);
}

PeakIonPairPtrVec findPairs(ExtendMsPtr ms,TheoPeakPtrVec ions,
                            int bgn,int end){
	PeakIonPairPtrVec pair_list;
	std::sort(ions.begin(),ions.end(),theo_peak_up);
	std::vector<double> ion_masses;
	for(unsigned int i=0;i<ions.size();i++){
		ion_masses.push_back(ions[i]->getModMass());
	}

	std::vector<double> real_masses;
	for(unsigned int i=0;i<ms->size();i++){
		real_masses.push_back(ms->getPeakPtr(i)->getPosition());
	}

	unsigned int i=0;
	unsigned int j=0;
	while(i<ms->size()&& j < ions.size()){
		double deviation = ms->getPeakPtr(i)->getPosition() 
        - ions[j]->getModMass();
		IonPtr ion = ions[j]->getIonPtr();
		double err = ms->getPeakPtr(i)->getOrigTolerance();
		if(ion->getPos()>=bgn && ion->getPos()<=end){
			if(std::abs(deviation)<=err){
				pair_list.push_back(
            PeakIonPairPtr(new PeakIonPair(ms->getPeakPtr(i),ions[j])));
			}
		}
		if(increaseIJ(i,j,deviation,err,real_masses,ion_masses)){
			i++;
		}
		else{
			j++;
		}
	}

	return pair_list;
}


std::vector<double> getNCScore(ExtendMsPtr ms,TheoPeakPtrVec ions,
                               int bgn,int end,double delta,double ppo){
	std::vector<double> ms_masses;
	for(unsigned int i=0;i<ms->size();i++){
		ms_masses.push_back(ms->getPeakPtr(i)->getPosition());
	}

	std::vector<double> theo_n_mass_list;

	for(unsigned int i=0;i<ions.size();i++){
		IonPtr ion = ions[i]->getIonPtr();
		int pos = ion->getPos();
		if(ion->getIonTypePtr()->isNTerm() && pos>=bgn && pos <= end){
			theo_n_mass_list.push_back(ions[i]->getModMass()+delta);
		}
	}

	double n_score = compUniqScore(ms_masses,theo_n_mass_list,ppo);

	std::vector<double> theo_c_mass_list;

	for(unsigned int i=0;i<ions.size();i++){
		IonPtr ion = ions[i]->getIonPtr();
		int pos = ion->getPos();
		if(!ion->getIonTypePtr()->isNTerm() && pos>=bgn && pos <= end){
			theo_c_mass_list.push_back(ions[i]->getModMass()+delta);
		}
	}

	double c_score =  compUniqScore(ms_masses,theo_c_mass_list,ppo);

	std::vector<double> result;
	result.push_back(n_score);
	result.push_back(c_score);
	return result;
}

} /* namespace prot */
