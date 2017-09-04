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

#include <numeric>
#include <algorithm>

#include "base/logger.hpp"
#include "spec/base_peak_type.hpp"
#include "spec/prm_ms_factory.hpp"

namespace prot {

void addTwoMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr, 
                  double prec_mono_mass, ActivationPtr active_type_ptr, PeakTolerancePtr tole_ptr) {
  double ori_mass = deconv_peak_ptr->getMonoMass();
  double n_term_mass = ori_mass - active_type_ptr->getNShift();
  PrmPeakPtr new_peak_ptr 
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,BasePeakType::ORIGINAL,n_term_mass,1);
  new_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  new_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  new_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  list.push_back(new_peak_ptr);
  double reverse_mass = prec_mono_mass - (deconv_peak_ptr->getMonoMass()-active_type_ptr->getCShift());
  PrmPeakPtr reverse_peak_ptr 
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,BasePeakType::REVERSED,reverse_mass,1);
  reverse_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  list.push_back(reverse_peak_ptr);
}

void addSuffixTwoMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr, 
                        double prec_mono_mass, ActivationPtr active_type_ptr, PeakTolerancePtr tole_ptr) {
  double ori_mass = deconv_peak_ptr->getMonoMass();
  double c_res_mass = ori_mass - active_type_ptr->getCShift() - MassConstant::getWaterMass();
  PrmPeakPtr new_peak_ptr 
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,BasePeakType::ORIGINAL,c_res_mass,1);
  new_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  new_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  new_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  list.push_back(new_peak_ptr);
  double reverse_mass = prec_mono_mass - (deconv_peak_ptr->getMonoMass()-active_type_ptr->getNShift())
      - MassConstant::getWaterMass();
  PrmPeakPtr reverse_peak_ptr 
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,BasePeakType::REVERSED,reverse_mass,1);
  reverse_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  list.push_back(reverse_peak_ptr);
}

void addSixMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr,
                  double prec_mono_mass, ActivationPtr active_type_ptr,
                  PeakTolerancePtr tole_ptr, const std::vector<double> &offsets){
  double ori_mass = deconv_peak_ptr->getMonoMass();
  for(size_t i = 0;i<offsets.size();i++){
    double mass = ori_mass - active_type_ptr->getNShift()+offsets[i];
    PrmPeakPtr peak_ptr = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,BasePeakType::ORIGINAL,mass,1);
    peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
    list.push_back(peak_ptr);
  }
  for(size_t i = 0;i<offsets.size();i++){
    double mass = prec_mono_mass-(ori_mass-active_type_ptr->getCShift()+offsets[i]);
    PrmPeakPtr peak_ptr = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,BasePeakType::REVERSED,mass,1);
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

PrmMsPtr geneMsTwoPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                      double prec_mono_mass){
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  double new_prec_mono_mass = prec_mono_mass - std::accumulate(sp_para_ptr->mod_mass_.begin(), sp_para_ptr->mod_mass_.end(), 0.0);
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, new_prec_mono_mass);
  //getSpTwoPrmPeak
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  PrmPeakPtrVec list;
  for(size_t i = 0; i < deconv_ms_ptr->size(); i++){
    addTwoMasses(list, spec_id,deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                 active_type_ptr, tole_ptr);
  }
  //filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());
  //sort 
  std::sort(list_filtered.begin(), list_filtered.end(), PrmPeak::cmpPosInc);
  int k = static_cast<int>(sp_para_ptr->mod_mass_.size()) + 1;
  for (size_t i = 0; i < list_filtered.size(); i++) {
    double m = list_filtered[i]->getMonoMass();
    for (int j = 1; j < k; j++) {
      if (list_filtered[i]->getMonoMass() > prec_mono_mass * j / k) {
        m -= sp_para_ptr->mod_mass_[j - 1];
      }
    }
    list_filtered[i] = std::make_shared<PrmPeak>(list_filtered[i]->getSpectrumId(),
                                                 list_filtered[i]->getBasePeakPtr(),
                                                 list_filtered[i]->getBaseTypePtr(),
                                                 m,
                                                 list_filtered[i]->getScore(),
                                                 list_filtered[i]->getStrictTolerance(),
                                                 list_filtered[i]->getNStrictCRelaxTolerance(),
                                                 list_filtered[i]->getNRelaxCStrictTolerance());
  }
  return PrmMsPtr(new Ms<PrmPeakPtr>(header_ptr,list_filtered)) ;
}

PrmMsPtr geneSuffixMsTwoPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                            double prec_mono_mass){
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  double new_prec_mono_mass = prec_mono_mass - std::accumulate(sp_para_ptr->mod_mass_.begin(), sp_para_ptr->mod_mass_.end(), 0.0);
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, new_prec_mono_mass);
  //getSpTwoPrmPeak
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  PrmPeakPtrVec list;
  for(size_t i = 0;i<deconv_ms_ptr->size();i++){
    addSuffixTwoMasses(list, spec_id,deconv_ms_ptr->getPeakPtr(i), new_prec_mono_mass,
                       active_type_ptr, tole_ptr);
  }
  //filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());
  //sort 
  std::sort(list_filtered.begin(), list_filtered.end(),PrmPeak::cmpPosInc);
  int k = static_cast<int>(sp_para_ptr->mod_mass_.size()) + 1;
  for (size_t i = 0; i < list_filtered.size(); i++) {
    double m = list_filtered[i]->getMonoMass();
    for (int j = 1; j < k; j++) {
      if (list_filtered[i]->getMonoMass() > prec_mono_mass * j / k) {
        m -= sp_para_ptr->mod_mass_[k - j -1]; 
      }
    }
    list_filtered[i] = std::make_shared<PrmPeak>(list_filtered[i]->getSpectrumId(),
                                                 list_filtered[i]->getBasePeakPtr(),
                                                 list_filtered[i]->getBaseTypePtr(),
                                                 m,
                                                 list_filtered[i]->getScore(),
                                                 list_filtered[i]->getStrictTolerance(),
                                                 list_filtered[i]->getNStrictCRelaxTolerance(),
                                                 list_filtered[i]->getNRelaxCStrictTolerance());
  }
  return PrmMsPtr(new Ms<PrmPeakPtr>(header_ptr,list_filtered)) ;
}

PrmMsPtr geneMsSixPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                      double prec_mono_mass){
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, prec_mono_mass);
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
    } else{
      addSixMasses(list,spec_id,deconv_ms_ptr->getPeakPtr(i),prec_mono_mass,
                   active_type_ptr, tole_ptr, sp_para_ptr->getExtendOffsets());
    }
  }

  //filterPrmPeak
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());
  //sort 
  std::sort(list_filtered.begin(), list_filtered.end(),PrmPeak::cmpPosInc);

  return PrmMsPtr(new Ms<PrmPeakPtr>(header_ptr,list_filtered)) ;
}

PrmMsPtr geneShiftMsSixPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr, 
                           double prec_mono_mass, double shift) {
  PrmMsPtr prm_ms_ptr = geneMsSixPtr(deconv_ms_ptr,spec_id, sp_para_ptr, prec_mono_mass);
  MsHeaderPtr ori_header_ptr = prm_ms_ptr->getMsHeaderPtr();
  MsHeaderPtr header_ptr = MsHeaderPtr(new MsHeader(*ori_header_ptr.get()));
  double mono_mz = (header_ptr->getPrecMonoMass()+shift) /header_ptr->getPrecCharge();
  header_ptr->setPrecMonoMz(mono_mz);
  PrmPeakPtrVec prm_peak_list ;
  for(size_t i=0;i< prm_ms_ptr->size();i++){
    double pos= prm_ms_ptr->getPeakPtr(i)->getPosition()+shift;
    if(pos>0){
      prm_ms_ptr->getPeakPtr(i)->setPosition(pos);
      prm_peak_list.push_back(prm_ms_ptr->getPeakPtr(i));
    }
  }
  return PrmMsPtr(new Ms<PrmPeakPtr>(header_ptr, prm_peak_list));
}

PrmMsPtrVec PrmMsFactory::geneMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                          SpParaPtr sp_para_ptr,
                                          double prec_mono_mass){
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneMsTwoPtr(deconv_ms_ptr_vec[i], i, 
                                          sp_para_ptr, prec_mono_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec PrmMsFactory::geneSuffixMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                                SpParaPtr sp_para_ptr,
                                                double prec_mono_mass){
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneSuffixMsTwoPtr(deconv_ms_ptr_vec[i], i, 
                                                sp_para_ptr, prec_mono_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec PrmMsFactory::geneMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                          SpParaPtr sp_para_ptr,
                                          double prec_mono_mass){
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneMsSixPtr(deconv_ms_ptr_vec[i], i,
                                          sp_para_ptr, prec_mono_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec PrmMsFactory::geneShiftMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                               SpParaPtr sp_para_ptr, double prec_mono_mass, 
                                               double shift) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneShiftMsSixPtr(deconv_ms_ptr_vec[i],i, sp_para_ptr,
                                               prec_mono_mass, shift));
  }
  return prm_ms_ptr_vec;
}

} /* namespace prot */
