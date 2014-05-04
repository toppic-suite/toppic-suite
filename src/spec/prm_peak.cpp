/*
 * prm_peak.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#include <algorithm>
#include <iostream>

#include "base/base_data.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

PrmPeak::PrmPeak(const DeconvPeakPtr &base_peak, int base_type,
                 double mono_mass, double score)
    :Peak(mono_mass,base_peak->getIntensity()){
      base_peak_=base_peak;
      base_type_=base_type;
      mono_mass_=mono_mass;
      score_=score;
      strict_tolerance_=0;
      n_strict_c_relax_tolerance_=0;
      n_relax_c_strict_tolerance_=0;
    }

void PrmPeak::addNghbEdge(const DeconvPeakPtr &peak,double offset,
                          const SPTypePtr &peak_type,double score){
  score_ +=score;
  SupportPeakPtr support_peak_ptr 
      = SupportPeakPtr(new SupportPeak(peak,offset,score,peak_type));
  neighbor_list_.push_back(support_peak_ptr);
}

int PrmPeak::getBreakType() {
  int break_type = 0;
  for(unsigned int i=0;i<neighbor_list_.size();i++){
    if(neighbor_list_[i]->getPeakTypePtr() == 
       SPTypeFactory::getSPTypePtr_N_TERM()){
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

/* 
 * in several local function, reference arguments are used as the 
 * output of function
 */
void addTwoMasses(PrmPeakPtrVec &list, const DeconvPeakPtr &peak,
                  double prec_mono_mass, ActivationPtr active_type){
  double orig_mass = peak->getMonoMass() - active_type->getNShift();
  PrmPeakPtr new_peak 
      = PrmPeakPtr(new PrmPeak(peak,PRM_PEAK_TYPE_ORIGINAL,orig_mass,1));
  list.push_back(new_peak);

  double reverse_mass = prec_mono_mass - 
      (peak->getMonoMass()-active_type->getCShift());
  PrmPeakPtr reverse_peak 
      = PrmPeakPtr(new PrmPeak(peak,PRM_PEAK_TYPE_REVERSED,reverse_mass,1));
  list.push_back(reverse_peak);
}

void addSixMasses(PrmPeakPtrVec &list, const DeconvPeakPtr &peak,
                  double prec_mono_mass, const ActivationPtr &active_type,
                  const std::vector<double> &offsets){
  for(unsigned int i = 0;i<offsets.size();i++){
    double mass = peak->getMonoMass()-active_type->getNShift()+offsets[i];
    list.push_back(PrmPeakPtr(
            new PrmPeak(peak,PRM_PEAK_TYPE_ORIGINAL,mass,1)
            ));
  }
  for(unsigned int i = 0;i<offsets.size();i++){
    double mass = prec_mono_mass-(peak->getMonoMass()-active_type->getCShift()+offsets[i]);
    list.push_back(PrmPeakPtr(
            new PrmPeak(peak,PRM_PEAK_TYPE_REVERSED,mass,1)
            ));
  }
}

void filterPeaks(const PrmPeakPtrVec &peak_list, PrmPeakPtrVec &filtered_list,
                 double prec_mono_mass, double min_mass) {
  for(unsigned int i =0;i<peak_list.size();i++){
    if(peak_list[i]->getPosition() >= min_mass &&
       peak_list[i]->getPosition()  <= prec_mono_mass - min_mass) {
      filtered_list.push_back(peak_list[i]);
    }
  }
}

void addZeroPrecPeak(PrmPeakPtrVec &peak_list, double prec_mono_mass) {
  //zero_peak
  DeconvPeakPtr zero_peak = DeconvPeakPtr(new DeconvPeak(-1,0,0,0));
  PrmPeakPtr prm_peak_ptr(new PrmPeak(zero_peak,PRM_PEAK_TYPE_ORIGINAL,0,1));
  peak_list.push_back(prm_peak_ptr);

  //prec_peak
  double prec_peak_shift = IonTypeFactory::getIonTypePtr_PREC()->getShift();
  double prec_peak_mass = prec_mono_mass - prec_peak_shift;
  DeconvPeakPtr prec_peak(new DeconvPeak(-1,prec_peak_mass, 0,0));
  prm_peak_ptr = PrmPeakPtr(
      new PrmPeak(prec_peak,PRM_PEAK_TYPE_ORIGINAL,prec_peak_mass, 1));
  peak_list.push_back(prm_peak_ptr);
  std::sort(peak_list.begin(), peak_list.end(),prmPeakUp);
}

void setTolerance(PrmPeakPtrVec &list, const PeakTolerancePtr &tole_ptr, 
                  double prec_mono_mass) {
  //peak with mass 0
  list[0]->setStrictTolerance(tole_ptr->compStrictErrorTole(0));
  list[0]->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(0));
  list[0]->setNRelaxCStrictTolerance(
      tole_ptr->compStrictErrorTole(prec_mono_mass));

  //peak with mass prec_mass -18
  int last_idx = list.size() -1;
  list[last_idx]->setStrictTolerance(
      tole_ptr->compStrictErrorTole(prec_mono_mass));
  list[last_idx]->setNStrictCRelacTolerance(
      tole_ptr->compStrictErrorTole(prec_mono_mass));
  list[last_idx]->setNRelaxCStrictTolerance(
      tole_ptr->compStrictErrorTole(0));
  //from 1 to size -1, not from 0 to size
  for(unsigned int i=1; i< list.size()-1;i++){
    double mass = list[i]->getBasePeak()->getMonoMass();
    list[i]->setStrictTolerance(tole_ptr->compStrictErrorTole(mass));
    if(list[i]->getBaseType()==PRM_PEAK_TYPE_ORIGINAL){
      list[i]->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(mass));
      list[i]->setNRelaxCStrictTolerance(
          tole_ptr->compRelaxErrorTole(mass, prec_mono_mass));
    }
    else{
      //if the peak is M-m ,use base peak and parent mass 
      //to compute the error tolerance
      list[i]->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(mass));
      list[i]->setNStrictCRelacTolerance(
          tole_ptr->compRelaxErrorTole(mass, prec_mono_mass));
    }
  }
}
    //end settolerance

PrmMsPtr getMsTwo(const DeconvMsPtr &deconv_ms,double delta, const SpParaPtr &sp_para){

  MsHeaderPtr header = prot::getDeltaHeaderPtr(deconv_ms,delta);

  //getSpTwoPrmPeak
  double prec_mono_mass = header->getPrecMonoMass();
  ActivationPtr active_type = header->getActivationPtr();
  PrmPeakPtrVec list;
  for(unsigned int i = 0;i<deconv_ms->size();i++){
    addTwoMasses(list,deconv_ms->getPeakPtr(i),prec_mono_mass,active_type);
  }
  //filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para->getMinMass());

  addZeroPrecPeak(list_filtered, prec_mono_mass);

  //setTolerance
  setTolerance(list_filtered, sp_para->getPeakTolerancePtr(), prec_mono_mass);

  double ppo = sp_para->getPeakTolerancePtr()->getPpo();
  return PrmMsPtr(new Ms<PrmPeakPtr>(header,list_filtered, ppo)) ;
}

PrmMsPtr getMsSix(const DeconvMsPtr &deconv_ms,double delta, 
                  const SpParaPtr &sp_para){

  MsHeaderPtr header = prot::getDeltaHeaderPtr(deconv_ms,delta);
  //getSpSixPrmPeak
  double prec_mono_mass = header->getPrecMonoMass();
  ActivationPtr active_type = header->getActivationPtr();
  double extend_min_mass = sp_para->getExtendMinMass();
  PrmPeakPtrVec list;
  //    std::cout<<deconv_ms->size()<<std::endl;
  for(unsigned int i=0;i< deconv_ms->size();i++){
    if(deconv_ms->getPeakPtr(i)->getMonoMass() <= extend_min_mass) {
      addTwoMasses(list,deconv_ms->getPeakPtr(i),prec_mono_mass,active_type);
    }
    else{
      addSixMasses(list,deconv_ms->getPeakPtr(i),prec_mono_mass,active_type,
                   sp_para->getExtendOffsets());
    }
  }


  //filterPrmPeak
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para->getMinMass());


  addZeroPrecPeak(list_filtered, prec_mono_mass);

  //setTolerance
  setTolerance(list_filtered, sp_para->getPeakTolerancePtr(), prec_mono_mass);

  //end settolerance
  double ppo = sp_para->getPeakTolerancePtr()->getPpo();
  return PrmMsPtr(new Ms<PrmPeakPtr>(header,list_filtered, ppo)) ;
}

PrmMsPtr getShiftMsSix(const DeconvMsPtr &deconv_ms,double delta,double shift,
                       const SpParaPtr &sp_para){
  PrmMsPtr ms = getMsSix(deconv_ms,delta,sp_para);
  MsHeaderPtr header_ptr = ms->getHeaderPtr();
  double mono_mz = (header_ptr->getPrecMonoMass()+shift)
      /header_ptr->getPrecCharge();
  header_ptr->setPrecMonoMz(mono_mz);
  PrmPeakPtrVec prm_peaks ;
  for(unsigned int i=0;i< ms->size();i++){
    double pos= ms->getPeakPtr(i)->getPosition()+shift;
    if(pos>0){
      ms->getPeakPtr(i)->setPosition(pos);
      prm_peaks.push_back(ms->getPeakPtr(i));
    }
  }
  double ppo = sp_para->getPeakTolerancePtr()->getPpo();
  return PrmMsPtr(new Ms<PrmPeakPtr>(ms->getHeaderPtr(),prm_peaks, ppo));
}

std::vector<std::vector<int>> getIntMassErrorList(const PrmMsPtr &ms,
                                                  double scale,
                                                  bool n_strict,
                                                  bool c_strict){
  std::vector<int> masses;
  std::vector<int> errors;
  int last_mass = -1;
  int last_error = 0;
  for(unsigned int i=0;i<ms->size();i++){
    int m = (int)std::round(ms->getPeakPtr(i)->getPosition()*scale);
    int e = 0;
    if(n_strict && c_strict){
      e = (int) std::ceil(ms->getPeakPtr(i)->getStrictTolerance()*scale);
    }
    else if(n_strict && !c_strict){
      e = (int) std::ceil(ms->getPeakPtr(i)->getNStrictCRelaxTolerance()*scale);
    }
    else if(!n_strict && c_strict){
      e = (int) std::ceil(ms->getPeakPtr(i)->getNRelaxCStrictTolerance()*scale);
    }
    if(m != last_mass){
      masses.push_back(m);
      errors.push_back(e);
      last_mass = m;
      last_error = e;
    }
    else if(e>last_error){
      errors.pop_back();
      errors.push_back(e);
      last_error = e;
    }
  }
  std::vector<int> mass_temp;
  std::vector<int> error_temp;
  std::vector<std::vector<int>> results;
  for(unsigned int i=0;i<masses.size();i++){
    mass_temp.push_back(masses[i]);
    error_temp.push_back(errors[i]);
  }
  results.push_back(mass_temp);
  results.push_back(error_temp);
  return results;
}

std::vector<double> getMassList(const PrmMsPtr &ms){
  std::vector<double> results;
  for(unsigned int i=0;i<ms->size();i++){
    results.push_back(ms->getPeakPtr(i)->getPosition());
  }
  return results;
}

std::vector<double> getScoreList(const PrmMsPtr &ms){
  std::vector<double> results;
  for(unsigned int i=0;i<ms->size();i++){
    results.push_back(ms->getPeakPtr(i)->getScr());
  }
  return results;
}

} /* namespace prot */
