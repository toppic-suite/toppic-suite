
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "base/logger.hpp"
#include "base/bp_spec.hpp"

namespace prot {

BpSpec::BpSpec(const ResSeqPtr &res_seq_ptr){
  initBreakPoints(res_seq_ptr);
}

void BpSpec::initBreakPoints(const ResSeqPtr &res_seq_ptr){
  int ext_len= res_seq_ptr->getLen()+1;
  if(ext_len <= 1){
    ext_len = 2;
  }
  // first breakpoint
  BreakPointPtr first_ptr(new BreakPoint(0,res_seq_ptr->getResMassSum()));
  break_point_ptr_vec_.push_back(first_ptr);

  double prm = 0;
  for(int i=0; i<res_seq_ptr->getLen()-1; i++){
    prm += res_seq_ptr->getResiduePtr(i)->getMass();
    double srm = res_seq_ptr->getResMassSum()-prm;
    if(srm <0){
      LOG_WARN("prms is larger than total mass! ");
    }
    break_point_ptr_vec_.push_back(BreakPointPtr(new BreakPoint(prm,srm)));
  }
  // last breakpoint
  BreakPointPtr last_ptr((new BreakPoint(res_seq_ptr->getResMassSum(),0)));
  break_point_ptr_vec_.push_back(last_ptr);
}

/* Get neutral ion masses for a specific ion type */
std::vector<double> BpSpec::getBreakPointMasses(IonTypePtr ion_type_ptr){
  std::vector<double> bp_mass_vec;
  if (ion_type_ptr->isNTerm()) {
    for (size_t i = 0; i < break_point_ptr_vec_.size(); i++) {
      bp_mass_vec.push_back(break_point_ptr_vec_[i]->getNTermMass(ion_type_ptr));
    }
  }
  else {
    for (size_t i = 0; i < break_point_ptr_vec_.size(); i++) {
      bp_mass_vec.push_back(break_point_ptr_vec_[i]->getCTermMass(ion_type_ptr));
    }
  }
  std::sort(bp_mass_vec.begin(),bp_mass_vec.end(),std::less<double>());
  return bp_mass_vec;
}

std::vector<double> BpSpec::getPrmMasses() {
  std::vector<double> mass_vec;
  for (size_t i = 0; i < break_point_ptr_vec_.size(); i++) {
    mass_vec.push_back(break_point_ptr_vec_[i]->getPrm());
  }
  std::sort(mass_vec.begin(),mass_vec.end(),std::less<double>());
  return mass_vec;
}

std::vector<double> BpSpec::getSrmMasses() {
  std::vector<double> mass_vec;
  for (size_t i = break_point_ptr_vec_.size() -1; i >= 0; i--) {
    mass_vec.push_back(break_point_ptr_vec_[i]->getSrm());
  }
  std::sort(mass_vec.begin(),mass_vec.end(),std::less<double>());
  return mass_vec;
}

/* Get rounded scaled neutral ion masses */ 
std::vector<int> BpSpec::getScaledMass(double scale, IonTypePtr ion_type_ptr){
  std::vector<int> result;
  if (ion_type_ptr->isNTerm()) {
    for(size_t i=0; i < break_point_ptr_vec_.size();i++){
      double value = break_point_ptr_vec_[i]->getNTermMass(ion_type_ptr)*scale;
      result.push_back(std::floor(value+0.5));
    }
  }
  else {
    for(size_t i=0; i < break_point_ptr_vec_.size();i++){
      double value = break_point_ptr_vec_[i]->getCTermMass(ion_type_ptr)*scale;
      result.push_back(std::floor(value+0.5));
    }
  }
  return result;
}

std::vector<int> BpSpec::getScaledPrmMasses(double scale){
  std::vector<int> result;
  for(size_t i=break_point_ptr_vec_.size()-1; i >= 0; i--){
    double value = break_point_ptr_vec_[i]->getPrm()*scale;
    result.push_back(std::floor(value+0.5));
  }
  return result;
}


std::vector<int> BpSpec::getScaledSrmMasses(double scale){
  std::vector<int> result;
  for(size_t i=break_point_ptr_vec_.size()-1; i >= 0; i--){
    double value = break_point_ptr_vec_[i]->getSrm()*scale;
    result.push_back(std::floor(value+0.5));
  }
  return result;
}

} /* namespace prot */
