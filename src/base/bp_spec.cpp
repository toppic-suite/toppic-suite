// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



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
  for (int i = break_point_ptr_vec_.size() -1; i >= 0; i--) {
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
  for(size_t i=0; i < break_point_ptr_vec_.size(); i++){
    double value = break_point_ptr_vec_[i]->getPrm()*scale;
    result.push_back(std::floor(value+0.5));
  }
  return result;
}


std::vector<int> BpSpec::getScaledSrmMasses(double scale){
  std::vector<int> result;
  for(int i=break_point_ptr_vec_.size()-1; i >= 0; i--){
    double value = break_point_ptr_vec_[i]->getSrm()*scale;
    result.push_back(std::floor(value+0.5));
  }
  return result;
}

} /* namespace prot */
