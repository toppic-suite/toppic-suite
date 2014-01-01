/*
 * pair.cpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

#include <prsm/pair.hpp>
#include <algorithm>
#include <float.h>

namespace prot {
bool increaseIJ(int i,int j,double deviation,double tolerace,std::vector<double> ms_masses,std::vector<double> seq_masses){
	if(deviation <=0){
		return true;
	}

	if(i>=ms_masses.size()-1){
		return false;
	}

	double next_pos = ms_masses[i+1];

	if(std::abs(next_pos-seq_masses[j])<= tolerace
			&& (j == seq_masses.size()-1||std::abs(next_pos - seq_masses[j]) < std::abs(next_pos - seq_masses[j+1]))){
		return true;
	}
	else{
		return false;
	}
}

std::vector<double> compPpoDeviation(std::vector<double> ms_masses,std::vector<double> theo_masses,double ppo){
	std::vector<double> min_distances;
	for(int i =0;i< ms_masses.size();i++){
		min_distances.push_back(DBL_MAX);
	}
	int i=0;
	int j=0;
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

	std::vector<double> min_ppos;

	for(int i=0;i<ms_masses.size();i++){
		min_ppos.push_back(DBL_MAX);
		if(ms_masses[i]>0){
			min_ppos[i] = min_distances[i]/ms_masses[i];
		}
	}
	return min_ppos;
}

double compUniqScore(std::vector<double> ms_masses,std::vector<double> theo_masses,double ppo){
	std::vector<double> min_distances;
	for(int i =0;i< ms_masses.size();i++){
		min_distances.push_back(DBL_MAX);
	}
	int i=0;
	int j=0;
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

	for(int i=0;i<theo_masses.size();i++){
		if(theo_masses[i]>0){
			double error = min_distances[i]/theo_masses[i];
			if(std::abs(error)<=ppo){
				score+=1.0;
			}
		}
	}

	return score;
}

} /* namespace prot */
