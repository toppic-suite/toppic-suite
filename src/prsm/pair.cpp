#include <limits>

#include "prsm/pair.hpp"

namespace prot {

Pair::Pair (int x, int y) {
  x_=x; y_=y;
}

std::vector<double> compPpoDeviation(const std::vector<double> &ms_masses,
                                     const std::vector<double> &theo_masses,
                                     double ppo){
  std::vector<double> min_distances;
  for(size_t i =0;i< ms_masses.size();i++){
    min_distances.push_back(std::numeric_limits<double>::infinity());
  }
  size_t i=0;
  size_t j=0;
  while(i < ms_masses.size() && j < theo_masses.size()) {
    double d = ms_masses[i] - theo_masses[j];
    if(std::abs(d) < std::abs(min_distances[i])){
      min_distances[i] = d;
    }
    double tolerance = ms_masses[i] * ppo;
    if(increaseIJ(i,j,d,tolerance,ms_masses,theo_masses)){
      i++;
    }
    else{
      j++;
    }
  }

  std::vector<double> min_ppos;

  for(size_t i=0; i<ms_masses.size(); i++){
    min_ppos.push_back(std::numeric_limits<double>::infinity());
    if(ms_masses[i]>0){
      min_ppos[i] = min_distances[i]/ms_masses[i];
    }
  }
  return min_ppos;
}

double compUniqScore(const std::vector<double> &ms_masses,
                     const std::vector<double> &theo_masses,
                     double ppo) {
  std::vector<double> min_distances;
  for(size_t i =0;i< ms_masses.size();i++){
    min_distances.push_back(std::numeric_limits<double>::infinity());
  }
  size_t i=0;
  size_t j=0;
  while(i<ms_masses.size() && j < theo_masses.size()){
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

  for(size_t i=0;i<theo_masses.size();i++){
    if(theo_masses[i]>0){
      double error = min_distances[i]/theo_masses[i];
      if(std::abs(error)<=ppo){
        score+=1.0;
      }
    }
  }

  return score;
}


std::vector<double> compPpoDeviation(ExtendMsPtr ms_ptr,
                                     const TheoPeakPtrVec &peak_ptrs, double ppo){
  std::vector<double> theo_masses;
  for(size_t i=0; i < peak_ptrs.size(); i++){
    theo_masses.push_back(peak_ptrs[i]->getModMass());
  }
  std::vector<double> ms_masses;
  for(size_t i=0; i<ms_ptr->size();i++){
    ms_masses.push_back(ms_ptr->getPeakPtr(i)->getPosition());
  }
  return compPpoDeviation(ms_masses, theo_masses, ppo);
}

double compIonScore(ExtendMsPtr ms_ptr, const TheoPeakPtrVec &peak_ptrs,
                    double recal,double ppo){
  std::vector<double> theo_masses;
  for(size_t i=0; i<peak_ptrs.size(); i++){
    theo_masses.push_back(peak_ptrs[i]->getModMass());
  }
  std::vector<double> ms_masses;
  for(size_t i=0;i<ms_ptr->size();i++){
    ms_masses.push_back(ms_ptr->getPeakPtr(i)->getPosition()*(1+recal));
  }
  return compUniqScore(ms_masses,theo_masses,ppo);
}

PeakIonPairPtrVec findPairs(ExtendMsPtr ms_ptr, TheoPeakPtrVec peak_ptrs,
                            int bgn, int end, double add_tolerance) {
  PeakIonPairPtrVec pair_list;
  std::sort(peak_ptrs.begin(), peak_ptrs.end(),theoPeakUp);
  std::vector<double> theo_masses;
  for(size_t i=0; i<peak_ptrs.size(); i++){
    theo_masses.push_back(peak_ptrs[i]->getModMass());
  }

  std::vector<double> ms_masses;
  for(size_t i=0;i<ms_ptr->size();i++){
    ms_masses.push_back(ms_ptr->getPeakPtr(i)->getPosition());
  }

  size_t i=0;
  size_t j=0;
  while(i<ms_ptr->size()&& j < peak_ptrs.size()){
    double deviation = ms_ptr->getPeakPtr(i)->getPosition() 
        - peak_ptrs[j]->getModMass();
    IonPtr ion = peak_ptrs[j]->getIonPtr();
    double err = ms_ptr->getPeakPtr(i)->getOrigTolerance() + add_tolerance;
    if(ion->getPos()>=bgn && ion->getPos()<=end){
      if(std::abs(deviation)<=err){
        pair_list.push_back(
            PeakIonPairPtr(new PeakIonPair(ms_ptr->getPeakPtr(i),peak_ptrs[j])));
      }
    }
    if(increaseIJ(i,j,deviation,err,ms_masses,theo_masses)){
      i++;
    }
    else{
      j++;
    }
  }

  return pair_list;
}

std::vector<double> getNCScore(ExtendMsPtr ms_ptr, const TheoPeakPtrVec &peak_ptrs,
                               int bgn,int end,double delta,double ppo){
  std::vector<double> ms_masses;
  for(size_t i=0; i<ms_ptr->size(); i++){
    ms_masses.push_back(ms_ptr->getPeakPtr(i)->getPosition());
  }

  std::vector<double> theo_n_masses;

  for(size_t i=0; i<peak_ptrs.size();i++){
    IonPtr ion = peak_ptrs[i]->getIonPtr();
    int pos = ion->getPos();
    if(ion->getIonTypePtr()->isNTerm() && pos>=bgn && pos <= end){
      theo_n_masses.push_back(peak_ptrs[i]->getModMass()+delta);
    }
  }

  double n_score = compUniqScore(ms_masses,theo_n_masses,ppo);

  std::vector<double> theo_c_masses;

  for(size_t i=0;i<peak_ptrs.size();i++){
    IonPtr ion = peak_ptrs[i]->getIonPtr();
    int pos = ion->getPos();
    if(!ion->getIonTypePtr()->isNTerm() && pos>=bgn && pos <= end){
      theo_c_masses.push_back(peak_ptrs[i]->getModMass()+delta);
    }
  }

  double c_score =  compUniqScore(ms_masses,theo_c_masses,ppo);

  std::vector<double> result;
  result.push_back(n_score);
  result.push_back(c_score);
  return result;
}

} /* namespace prot */
