#include <algorithm>
#include <iostream>

#include "base/base_data.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

PrmPeak::PrmPeak(DeconvPeakPtr base_peak_ptr, int base_type,
                 double mono_mass, double score)
    :Peak(mono_mass, base_peak_ptr->getIntensity()){
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

/* 
 *  * in several local function, reference arguments are used as the 
 *   * output of function
 *    */
void addTwoMasses(PrmPeakPtrVec &list, DeconvPeakPtr peak_ptr,
                  double prec_mono_mass, ActivationPtr active_type_ptr){
  double orig_mass = peak_ptr->getMonoMass() - active_type_ptr->getNShift();
  PrmPeakPtr new_peak_ptr 
      = PrmPeakPtr(new PrmPeak(peak_ptr,PRM_PEAK_TYPE_ORIGINAL,orig_mass,1));
  list.push_back(new_peak_ptr);

  double reverse_mass = prec_mono_mass - 
      (peak_ptr->getMonoMass()-active_type_ptr->getCShift());
  PrmPeakPtr reverse_peak_ptr 
      = PrmPeakPtr(new PrmPeak(peak_ptr,PRM_PEAK_TYPE_REVERSED,reverse_mass,1));
  list.push_back(reverse_peak_ptr);
}

void addSixMasses(PrmPeakPtrVec &list, DeconvPeakPtr peak_ptr,
                  double prec_mono_mass, ActivationPtr active_type_ptr,
                  const std::vector<double> &offsets){
  for(size_t i = 0;i<offsets.size();i++){
    double mass = peak_ptr->getMonoMass()-active_type_ptr->getNShift()+offsets[i];
    list.push_back(PrmPeakPtr(
            new PrmPeak(peak_ptr,PRM_PEAK_TYPE_ORIGINAL,mass,1)
            ));
  }
  for(size_t i = 0;i<offsets.size();i++){
    double mass = prec_mono_mass-(peak_ptr->getMonoMass()-active_type_ptr->getCShift()+offsets[i]);
    list.push_back(PrmPeakPtr(
            new PrmPeak(peak_ptr,PRM_PEAK_TYPE_REVERSED,mass,1)
            ));
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

void addZeroPrecPeak(PrmPeakPtrVec &peak_list, double prec_mono_mass) {
  //zero_peak
  DeconvPeakPtr zero_peak_ptr = DeconvPeakPtr(new DeconvPeak(-1,0,0,0));
  PrmPeakPtr prm_peak_ptr(new PrmPeak(zero_peak_ptr,PRM_PEAK_TYPE_ORIGINAL,0,1));
  peak_list.push_back(prm_peak_ptr);

  //prec_peak
  double prec_peak_shift = IonTypeFactory::getIonTypePtr_PREC()->getShift();
  double prec_peak_mass = prec_mono_mass - prec_peak_shift;
  DeconvPeakPtr prec_peak_ptr(new DeconvPeak(-1,prec_peak_mass, 0,0));
  prm_peak_ptr = PrmPeakPtr(
      new PrmPeak(prec_peak_ptr,PRM_PEAK_TYPE_ORIGINAL,prec_peak_mass, 1));
  peak_list.push_back(prm_peak_ptr);
  std::sort(peak_list.begin(), peak_list.end(),prmPeakUp);
}

void setTolerance(PrmPeakPtrVec &list, PeakTolerancePtr tole_ptr, 
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
  for(size_t i=1; i< list.size()-1;i++){
    double mass = list[i]->getBasePeakPtr()->getMonoMass();
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

PrmMsPtrVec createMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                              SpParaPtr sp_para_ptr,
                              double prec_mono_mass){
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(createMsTwoPtr(deconv_ms_ptr_vec[i], sp_para_ptr,
                                            prec_mono_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec createMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                              SpParaPtr sp_para_ptr,
                              double prec_mono_mass){
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(createMsSixPtr(deconv_ms_ptr_vec[i], sp_para_ptr,
                                            prec_mono_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec createShiftMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                   SpParaPtr sp_para_ptr, double prec_mono_mass, 
                                   double shift) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(createShiftMsSixPtr(deconv_ms_ptr_vec[i], sp_para_ptr,
                                                 prec_mono_mass, shift));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtr createMsTwoPtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr,
                        double prec_mono_mass){

  MsHeaderPtr header_ptr = getHeaderPtr(deconv_ms_ptr, prec_mono_mass);

  //getSpTwoPrmPeak
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PrmPeakPtrVec list;
  for(size_t i = 0;i<deconv_ms_ptr->size();i++){
    addTwoMasses(list,deconv_ms_ptr->getPeakPtr(i),prec_mono_mass,active_type_ptr);
  }
  //filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());

  addZeroPrecPeak(list_filtered, prec_mono_mass);

  //setTolerance
  setTolerance(list_filtered, sp_para_ptr->getPeakTolerancePtr(), prec_mono_mass);

  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  return PrmMsPtr(new Ms<PrmPeakPtr>(header_ptr,list_filtered, ppo)) ;
}

PrmMsPtr createMsSixPtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr,
                        double prec_mono_mass){

  MsHeaderPtr header_ptr = getHeaderPtr(deconv_ms_ptr, prec_mono_mass);
  //getSpSixPrmPeak
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  double extend_min_mass = sp_para_ptr->getExtendMinMass();
  PrmPeakPtrVec list;
  //    std::cout<<deconv_ms->size()<<std::endl;
  for(size_t i=0;i< deconv_ms_ptr->size();i++){
    if(deconv_ms_ptr->getPeakPtr(i)->getMonoMass() <= extend_min_mass) {
      addTwoMasses(list,deconv_ms_ptr->getPeakPtr(i),prec_mono_mass,active_type_ptr);
    }
    else{
      addSixMasses(list,deconv_ms_ptr->getPeakPtr(i),prec_mono_mass,active_type_ptr,
                   sp_para_ptr->getExtendOffsets());
    }
  }


  //filterPrmPeak
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());


  addZeroPrecPeak(list_filtered, prec_mono_mass);

  //setTolerance
  setTolerance(list_filtered, sp_para_ptr->getPeakTolerancePtr(), prec_mono_mass);

  //end settolerance
  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  return PrmMsPtr(new Ms<PrmPeakPtr>(header_ptr,list_filtered, ppo)) ;
}

PrmMsPtr createShiftMsSixPtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr, 
                             double prec_mono_mass, double shift) {
  PrmMsPtr prm_ms_ptr = createMsSixPtr(deconv_ms_ptr,sp_para_ptr, prec_mono_mass);
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

inline bool massErrorUp(const std::pair<int, int> &a, const std::pair<int,int> b) {
  return a.first < b.first;
}

std::vector<std::pair<int, int>> getIntMassErrorList(const PrmMsPtrVec &prm_ms_ptr_vec, 
                                                     double scale, bool n_strict, bool c_strict){
  std::vector<std::pair<int,int>> mass_errors;
  for (size_t i = 0; i < prm_ms_ptr_vec.size(); i++) {
    PrmMsPtr prm_ms_ptr = prm_ms_ptr_vec[i];
    std::pair<int,int> last_mass_error(-1, 0);
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
      if(m != last_mass_error.first){
        mass_errors.push_back(std::pair<int,int>(m,e));
        last_mass_error = std::pair<int,int>(m,e);
      }
      else if(e>last_mass_error.second){
        mass_errors.pop_back();
        mass_errors.push_back(std::pair<int,int>(m,e));
        last_mass_error = std::pair<int,int>(m,e);
      }
    }
  }
  std::sort(mass_errors.begin(), mass_errors.end(),massErrorUp);
  return mass_errors;
}

std::vector<double> getMassList(PrmMsPtr prm_ms_ptr){
  std::vector<double> results;
  for(size_t i=0;i<prm_ms_ptr->size();i++){
    results.push_back(prm_ms_ptr->getPeakPtr(i)->getPosition());
  }
  return results;
}

std::vector<double> getScoreList(PrmMsPtr prm_ms_ptr){
  std::vector<double> results;
  for(size_t i=0;i<prm_ms_ptr->size();i++){
    results.push_back(prm_ms_ptr->getPeakPtr(i)->getScore());
  }
  return results;
}

} /* namespace prot */
