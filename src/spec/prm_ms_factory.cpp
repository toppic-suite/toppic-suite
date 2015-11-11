#include <algorithm>
#include <iostream>

#include "base/logger.hpp"
#include "base/base_data.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

PrmPeak::PrmPeak(int spec_id, DeconvPeakPtr base_peak_ptr, int base_type,
                 double mono_mass, double score)
    :Peak(mono_mass, base_peak_ptr->getIntensity()){
      spec_id_ = spec_id;
      base_peak_ptr_=base_peak_ptr;
      base_type_=base_type;
      mono_mass_=mono_mass;
      score_=score;
      strict_tolerance_=0;
      n_strict_c_relax_tolerance_=0;
      n_relax_c_strict_tolerance_=0;
    }

void PrmPeak::addNghbEdge(DeconvPeakPtr deconv_peak_ptr,double offset,
                          SPTypePtr peak_type_ptr,double score){
  score_ +=score;
  SupportPeakPtr support_peak_ptr 
      = SupportPeakPtr(new SupportPeak(deconv_peak_ptr,offset,score,peak_type_ptr));
  neighbor_list_.push_back(support_peak_ptr);
}

int PrmPeak::getBreakType() {
  int break_type = 0;
  for(size_t i=0; i<neighbor_list_.size(); i++){
    if(neighbor_list_[i]->getPeakTypePtr() == 
       SPTypeFactory::getSPTypePtr_N_TERM()) {
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

void addTwoMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr, 
                  double prec_mono_mass, ActivationPtr active_type_ptr, PeakTolerancePtr tole_ptr) {
  double ori_mass = deconv_peak_ptr->getMonoMass();
  double n_term_mass = ori_mass - active_type_ptr->getNShift();
  PrmPeakPtr new_peak_ptr 
      = PrmPeakPtr(new PrmPeak(spec_id, deconv_peak_ptr,PRM_PEAK_TYPE_ORIGINAL,n_term_mass,1));
  new_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  new_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  new_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  list.push_back(new_peak_ptr);
  double reverse_mass = prec_mono_mass - (deconv_peak_ptr->getMonoMass()-active_type_ptr->getCShift());
  PrmPeakPtr reverse_peak_ptr 
      = PrmPeakPtr(new PrmPeak(spec_id, deconv_peak_ptr,PRM_PEAK_TYPE_REVERSED,reverse_mass,1));
  reverse_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  list.push_back(reverse_peak_ptr);
}

void addSixMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr,
                  double prec_mono_mass, ActivationPtr active_type_ptr,
                  PeakTolerancePtr tole_ptr, const std::vector<double> &offsets){
  double ori_mass = deconv_peak_ptr->getMonoMass();
  for(size_t i = 0;i<offsets.size();i++){
    double mass = ori_mass - active_type_ptr->getNShift()+offsets[i];
    PrmPeakPtr peak_ptr = PrmPeakPtr(new PrmPeak(spec_id, deconv_peak_ptr,PRM_PEAK_TYPE_ORIGINAL,mass,1));
    peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
    list.push_back(peak_ptr);
  }
  for(size_t i = 0;i<offsets.size();i++){
    double mass = prec_mono_mass-(ori_mass-active_type_ptr->getCShift()+offsets[i]);
    PrmPeakPtr peak_ptr = PrmPeakPtr(new PrmPeak(spec_id, deconv_peak_ptr,PRM_PEAK_TYPE_REVERSED,mass,1));
    peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
    list.push_back(peak_ptr);
  }
}

void filterPeaks(const PrmPeakPtrVec &peak_list, PrmPeakPtrVec &filtered_list,
                 double prec_mono_mass, double min_mass) {
  for(size_t i =0;i<peak_list.size();i++){
    if(peak_list[i]->getPosition() >= min_mass &&
       peak_list[i]->getPosition()  <= prec_mono_mass - min_mass) {
      filtered_list.push_back(peak_list[i]);
    }
  }
}

PrmPeakPtr getZeroPeakPtr(int spec_id, double prec_mono_mass, PeakTolerancePtr tole_ptr, double score) {
  //zero_peak
  DeconvPeakPtr zero_peak_ptr = DeconvPeakPtr(new DeconvPeak(-1,0,0,0));
  PrmPeakPtr prm_peak_ptr(new PrmPeak(spec_id, zero_peak_ptr,PRM_PEAK_TYPE_ORIGINAL,0, score));
  // set tolerance 
  prm_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(0));
  prm_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(0));
  prm_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  return prm_peak_ptr;
}

PrmPeakPtr getPrecPeakPtr(int spec_id, double prec_mono_mass, PeakTolerancePtr tole_ptr, double score) {
  //prec_peak
  double prec_peak_shift = IonTypeFactory::getIonTypePtr_PREC()->getShift();
  double prec_peak_mass = prec_mono_mass - prec_peak_shift;
  DeconvPeakPtr prec_peak_ptr(new DeconvPeak(-1,prec_peak_mass, 0,0));
  PrmPeakPtr prm_peak_ptr = PrmPeakPtr(
      new PrmPeak(spec_id,prec_peak_ptr,PRM_PEAK_TYPE_ORIGINAL,prec_peak_mass, score));
  // set tolerance
  prm_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  prm_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  prm_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(0));
  return prm_peak_ptr;
}

PrmMsPtr createMsTwoPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                        double prec_mono_mass){
  MsHeaderPtr header_ptr = getHeaderPtr(deconv_ms_ptr, prec_mono_mass);
  //getSpTwoPrmPeak
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  PrmPeakPtrVec list;
  for(size_t i = 0;i<deconv_ms_ptr->size();i++){
    addTwoMasses(list, spec_id,deconv_ms_ptr->getPeakPtr(i),prec_mono_mass,
                 active_type_ptr, tole_ptr);
  }
  //filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());
  //sort 
  std::sort(list_filtered.begin(), list_filtered.end(),prmPeakUp);

  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  return PrmMsPtr(new Ms<PrmPeakPtr>(header_ptr,list_filtered, ppo)) ;
}

PrmMsPtr createMsSixPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                        double prec_mono_mass){
  MsHeaderPtr header_ptr = getHeaderPtr(deconv_ms_ptr, prec_mono_mass);
  //getSpSixPrmPeak
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  double extend_min_mass = sp_para_ptr->getExtendMinMass();
  PrmPeakPtrVec list;
  //    std::cout<<deconv_ms->size()<<std::endl;
  for(size_t i=0;i< deconv_ms_ptr->size();i++){
    if(deconv_ms_ptr->getPeakPtr(i)->getMonoMass() <= extend_min_mass) {
      addTwoMasses(list, spec_id,deconv_ms_ptr->getPeakPtr(i),prec_mono_mass,
                   active_type_ptr, tole_ptr);
    }
    else{
      addSixMasses(list,spec_id,deconv_ms_ptr->getPeakPtr(i),prec_mono_mass,
                   active_type_ptr, tole_ptr, sp_para_ptr->getExtendOffsets());
    }
  }

  //filterPrmPeak
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());
  //sort 
  std::sort(list_filtered.begin(), list_filtered.end(),prmPeakUp);

  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  return PrmMsPtr(new Ms<PrmPeakPtr>(header_ptr,list_filtered, ppo)) ;
}

PrmMsPtr createShiftMsSixPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr, 
                             double prec_mono_mass, double shift) {
  PrmMsPtr prm_ms_ptr = createMsSixPtr(deconv_ms_ptr,spec_id, sp_para_ptr, prec_mono_mass);
  MsHeaderPtr ori_header_ptr = prm_ms_ptr->getHeaderPtr();
  MsHeaderPtr header_ptr = MsHeaderPtr(new MsHeader(*ori_header_ptr.get()));
  double mono_mz = (header_ptr->getPrecMonoMass()+shift)
      /header_ptr->getPrecCharge();
  header_ptr->setPrecMonoMz(mono_mz);
  PrmPeakPtrVec prm_peak_list ;
  for(size_t i=0;i< prm_ms_ptr->size();i++){
    double pos= prm_ms_ptr->getPeakPtr(i)->getPosition()+shift;
    if(pos>0){
      prm_ms_ptr->getPeakPtr(i)->setPosition(pos);
      prm_peak_list.push_back(prm_ms_ptr->getPeakPtr(i));
    }
  }
  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  return PrmMsPtr(new Ms<PrmPeakPtr>(header_ptr, prm_peak_list, ppo));
}

PrmMsPtrVec createMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                              SpParaPtr sp_para_ptr,
                              double prec_mono_mass){
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(createMsTwoPtr(deconv_ms_ptr_vec[i], i, 
                                            sp_para_ptr, prec_mono_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec createMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                              SpParaPtr sp_para_ptr,
                              double prec_mono_mass){
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(createMsSixPtr(deconv_ms_ptr_vec[i], i,
                                            sp_para_ptr, prec_mono_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec createShiftMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                   SpParaPtr sp_para_ptr, double prec_mono_mass, 
                                   double shift) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(createShiftMsSixPtr(deconv_ms_ptr_vec[i],i, sp_para_ptr,
                                                 prec_mono_mass, shift));
  }
  return prm_ms_ptr_vec;
}

inline bool massErrorUp(const std::pair<int, int> &a, const std::pair<int,int> b) {
  return a.first < b.first;
}

std::pair<int, int> getMassError(PrmPeakPtr peak_ptr, double scale, bool n_strict, bool c_strict) {
  int m = (int)std::round(peak_ptr->getPosition()*scale);
  int e = 0;
  if(n_strict && c_strict){
    e = (int) std::ceil(peak_ptr->getStrictTolerance()*scale);
  }
  else if(n_strict && !c_strict){
    e = (int) std::ceil(peak_ptr->getNStrictCRelaxTolerance()*scale);
  }
  else if(!n_strict && c_strict){
    e = (int) std::ceil(peak_ptr->getNRelaxCStrictTolerance()*scale);
  }
  std::pair<int, int> mass_error(m, e);
  return mass_error;
}

std::vector<std::pair<int, int>> getIntMassErrorList(const PrmMsPtrVec &prm_ms_ptr_vec, 
                                                     PeakTolerancePtr tole_ptr,
                                                     double scale, bool n_strict, bool c_strict){
  std::vector<std::pair<int,int>> mass_errors;
  for (size_t i = 0; i < prm_ms_ptr_vec.size(); i++) {
    PrmMsPtr prm_ms_ptr = prm_ms_ptr_vec[i];
    std::pair<int,int> last_mass_error(-1, 0);
    for(size_t j=0; j<prm_ms_ptr->size(); j++){
      std::pair<int, int> cur_m_e = getMassError(prm_ms_ptr->getPeakPtr(j), scale, n_strict, c_strict);
      if(cur_m_e.first != last_mass_error.first){
        mass_errors.push_back(cur_m_e);
        last_mass_error = cur_m_e;
      }
      else if(cur_m_e.second > last_mass_error.second){
        mass_errors.pop_back();
        mass_errors.push_back(cur_m_e);
        last_mass_error = cur_m_e;
      }
    }
    //add zero mass for each spectrum to increase the score for zero mass
    double prec_mass = prm_ms_ptr_vec[i]->getHeaderPtr()->getPrecMonoMass();
    PrmPeakPtr zero_prm_ptr = getZeroPeakPtr(i, prec_mass, tole_ptr, 1);
    mass_errors.push_back(getMassError(zero_prm_ptr, scale, n_strict, c_strict));
    //add prec mass for each spectrum 
    PrmPeakPtr prec_prm_ptr = getPrecPeakPtr(i, prec_mass, tole_ptr, 1);
    mass_errors.push_back(getMassError(prec_prm_ptr, scale, true, true));
  }

  std::sort(mass_errors.begin(), mass_errors.end(),massErrorUp);
  return mass_errors;
}

PrmPeakPtrVec getPrmPeakPtrs(const PrmMsPtrVec &prm_ms_ptr_vec, PeakTolerancePtr tole_ptr) {
  PrmPeakPtrVec peak_list;
  for (size_t i = 0; i < prm_ms_ptr_vec.size(); i++) {
    for(size_t j= 0;j<prm_ms_ptr_vec[i]->size() ;j++){
      peak_list.push_back(prm_ms_ptr_vec[i]->getPeakPtr(j));
    }
  }
  // add zero 
  double prec_mass = prm_ms_ptr_vec[0]->getHeaderPtr()->getPrecMonoMass();
  // use spec_id = 0 and score = group_spec_num (size of prm_ms_ptr_vec)
  PrmPeakPtr zero_prm_ptr = getZeroPeakPtr(0, prec_mass, tole_ptr, prm_ms_ptr_vec.size());
  peak_list.push_back(zero_prm_ptr);
  // add prec mass  
  PrmPeakPtr prec_prm_ptr = getPrecPeakPtr(0, prec_mass, tole_ptr, prm_ms_ptr_vec.size());
  peak_list.push_back(prec_prm_ptr);
  std::sort(peak_list.begin(), peak_list.end(), prmPeakUp);
  for (size_t i = 0; i < peak_list.size(); i++) {
    peak_list[i]->setPeakId(i);
  }
  return peak_list;
}


/*
void setTolerance(PrmPeakPtrVec &list, PeakTolerancePtr tole_ptr, 
                  double prec_mono_mass) {
  for(size_t i= 0; i< list.size();i++){
    double mass = list[i]->getBasePeakPtr()->getMonoMass();
    list[i]->setStrictTolerance(tole_ptr->compStrictErrorTole(mass));
    if(list[i]->getBaseType()==PRM_PEAK_TYPE_ORIGINAL){
      list[i]->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(mass));
      list[i]->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(mass, prec_mono_mass));
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
*/





/*
std::pair<std::vector<int>, std::vector<int>> getIntMassErrorList(PrmMsPtr prm_ms_ptr, double scale,
                                                                  bool n_strict, bool c_strict){
  std::vector<int> masses;
  std::vector<int> errors;
  int last_mass = -1;
  int last_error = 0;
  for(size_t i=0;i<prm_ms_ptr->size();i++){
    int m = (int)std::round(prm_ms_ptr->getPeakPtr(i)->getPosition()*scale);
    int e = 0;
    if(n_strict && c_strict){
      e = (int) std::ceil(prm_ms_ptr->getPeakPtr(i)->getStrictTolerance()*scale);
    }
    else if(n_strict && !c_strict){
      e = (int) std::ceil(prm_ms_ptr->getPeakPtr(i)->getNStrictCRelaxTolerance()*scale);
    }
    else if(!n_strict && c_strict){
      e = (int) std::ceil(prm_ms_ptr->getPeakPtr(i)->getNRelaxCStrictTolerance()*scale);
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
  std::pair<std::vector<int>, std::vector<int>> results( masses, errors);
  return results;
}
*/


/*
std::vector<double> getMassList(PrmMsPtr prm_ms_ptr){
  std::vector<double> results;
  for(size_t i=0;i<prm_ms_ptr->size();i++){
    results.push_back(prm_ms_ptr->getPeakPtr(i)->getPosition());
  }
  return results;
}

std::vector<double> getMassList(const PrmMsPtrVec &prm_ms_ptr_vec){
  std::vector<double> results;
  //add 0
  results.push_back(0);
  for (size_t i = 0; i < prm_ms_ptr_vec.size(); i++) {
    // skip mass 0 and prec_mass;
    for(size_t j=1;j<prm_ms_ptr_vec[i]->size() - 1 ;j++){
      results.push_back(prm_ms_ptr_vec[i]->getPeakPtr(j)->getPosition());
    }
  }
  // add prec mass
  double prec_mass = prm_ms_ptr_vec[0]->getPeakPtr(prm_ms_ptr_vec[0]->size()-1)->getPosition();
  results.push_back(prec_mass);
  std::sort(results.begin(), results.end(), std::less<double>());
  return results;
}

std::vector<double> getScoreList(PrmMsPtr prm_ms_ptr){
  std::vector<double> results;
  for(size_t i=0;i<prm_ms_ptr->size();i++){
    results.push_back(prm_ms_ptr->getPeakPtr(i)->getScore());
  }
  return results;
}
*/

} /* namespace prot */
