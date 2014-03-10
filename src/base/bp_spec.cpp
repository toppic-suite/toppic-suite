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

#include "base/logger.hpp"
#include "base/bp_spec.hpp"

namespace prot {

BpSpec::BpSpec(ResSeqPtr res_seq_ptr){
  initBreakPoints(res_seq_ptr);
}

void BpSpec::initBreakPoints(ResSeqPtr res_seq_ptr){
  seq_mass_ = res_seq_ptr->getSeqMass();
  int ext_len= res_seq_ptr->getLen()+1;
  if(ext_len <= 1){
    ext_len = 2;
  }
  BreakPointPtr first_ptr(new BreakPoint(0,res_seq_ptr->getResMassSum()));
  break_point_ptr_vec_.push_back(first_ptr);

  double prm = 0;
  for(int i=0; i<res_seq_ptr->getLen()-1; i++){
    prm += res_seq_ptr->getResiduePtr(i)->getMass();
    double srm = res_seq_ptr->getResMassSum()-prm;
    if(srm <0){
      LOG_WARN("prms is larger than totle mass! ");
    }
    break_point_ptr_vec_.push_back(BreakPointPtr(new BreakPoint(prm,srm)));
  }
  BreakPointPtr last_ptr((new BreakPoint(res_seq_ptr->getResMassSum(),0)));
  break_point_ptr_vec_.push_back(last_ptr);
}

//
std::vector<double> BpSpec::getBreakPointMasses(IonTypePtr ion_type_ptr){
  std::vector<double> bp_mass_vec;
  if (ion_type_ptr->isNTerm()) {
    for (unsigned int i = 0; i < break_point_ptr_vec_.size(); i++) {
      bp_mass_vec.push_back(break_point_ptr_vec_[i]->getNTermMass(ion_type_ptr));
    }
  }
  else {
    for (unsigned int i = 0; i < break_point_ptr_vec_.size(); i++) {
      bp_mass_vec.push_back(break_point_ptr_vec_[i]->getCTermMass(ion_type_ptr));
    }
  }
  std::sort(bp_mass_vec.begin(),bp_mass_vec.end(),std::less<double>());
  return bp_mass_vec;
}

//
std::vector<double> BpSpec::getPrmMasses() {
  std::vector<double> mass_vec;
  for (unsigned int i = 0; i < break_point_ptr_vec_.size(); i++) {
    mass_vec.push_back(break_point_ptr_vec_[i]->getPrm());
  }
  std::sort(mass_vec.begin(),mass_vec.end(),std::less<double>());
  return mass_vec;
}


void BpSpec::addBreakPointMass(double mass,double seq_mass,double min_mass,
                               std::vector<double> &mass_vec){
  if (mass >= min_mass && mass <= seq_mass - min_mass){
    mass_vec.push_back(mass);
  }
}

std::vector<double> BpSpec::getBreakPointMasses(double n_term_shift,
                                                double c_term_shift,
                                                double min_mass,
                                                IonTypePtr n_ion_type_ptr,
                                                IonTypePtr c_ion_type_ptr){
  std::vector<double> result;
  result.push_back(0.0);
  double new_seq_mass = seq_mass_ + n_term_shift + c_term_shift;
  //n
  for(unsigned int i=0; i < break_point_ptr_vec_.size();i++){
    double mass = break_point_ptr_vec_[i]->getNTermMass(n_ion_type_ptr)+n_term_shift;
    addBreakPointMass(mass,new_seq_mass,min_mass,result);
  }
  //c
  for(unsigned int i=0; i < break_point_ptr_vec_.size();i++){
    double mass = break_point_ptr_vec_[i]->getCTermMass(c_ion_type_ptr)+n_term_shift;
    addBreakPointMass(mass,new_seq_mass,min_mass,result);
  }
  result.push_back(new_seq_mass);
  std::sort(result.begin(),result.end(),std::less<double>());
  return result;
}

std::vector<int> BpSpec::getScaledMass(double scale,IonTypePtr ion_type_ptr){
  std::vector<int> result;
  if (ion_type_ptr->isNTerm()) {
    for(unsigned int i=0; i < break_point_ptr_vec_.size();i++){
      double value = break_point_ptr_vec_[i]->getNTermMass(ion_type_ptr)*scale;
      result.push_back(std::floor(value+0.5));
    }
  }
  else {
    for(unsigned int i=0; i < break_point_ptr_vec_.size();i++){
      double value = break_point_ptr_vec_[i]->getCTermMass(ion_type_ptr)*scale;
      result.push_back(std::floor(value+0.5));
    }
  }
  return result;
}

std::vector<int> BpSpec::getScaledPrmMasses(double scale){
  std::vector<int> result;
  for(unsigned int i=0; i < break_point_ptr_vec_.size();i++){
    double value = break_point_ptr_vec_[i]->getPrm()*scale;
    result.push_back(std::floor(value+0.5));
  }
  return result;
}

void BpSpec::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("bp_spec");
  std::string str = convertToString(seq_mass_);
  xml_doc->addElement(element, "seq_mass", str.c_str());
  xercesc::DOMElement* bplist = xml_doc->createElement("break_point_list");
  for(unsigned int i=0;i<break_point_ptr_vec_.size();i++){
    break_point_ptr_vec_[i]->appendXml(xml_doc,bplist);
  }
  element->appendChild(bplist);
  parent->appendChild(element);
}

int getFirstResPos(double n_term_shift,std::vector<double> ext_b_masses){
  double trunc_len = - n_term_shift;
  int best_pos = -1;
  double best_shift = std::numeric_limits<double>::infinity();
  for(unsigned int i = 0; i < ext_b_masses.size();i++){
    if(std::abs(ext_b_masses[i] - trunc_len) < best_shift){
      best_pos = i;
      best_shift = std::abs(ext_b_masses[i] - trunc_len);
    }
  }
  return best_pos;
}

int getLastResPos(double c_term_shift,std::vector<double> ext_b_masses){
  double trunc_len = -c_term_shift;
  int best_pos = -1;
  double best_shift = std::numeric_limits<double>::infinity();
  double pep_mass = ext_b_masses[ext_b_masses.size()-1];
  for(unsigned int i=0;i<ext_b_masses.size();i++){
    if (std::abs(pep_mass-ext_b_masses[i]-trunc_len)<best_shift){
      best_pos=i;
      best_shift = std::abs(pep_mass-ext_b_masses[i]-trunc_len);
    }
  }
  if(best_pos < 0){
    LOG_ERROR("get last residue position error! ");
    throw "get last residue position error!";
  }
  return best_pos - 1;
}

} /* namespace prot */
