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


#include "spec/extend_ms_factory.hpp"

namespace prot {

ExtendMsPtr ExtendMsFactory::geneMsThreePtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr, 
                                            double new_prec_mass) {
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, new_prec_mass);

  ExtendPeakPtrVec list;
  double ext_min_mass = sp_para_ptr->getExtendMinMass();
  std::vector<double> ext_offsets = sp_para_ptr->getExtendOffsets();
  int k = static_cast<int>(sp_para_ptr->mod_mass_.size()) + 1;
  for(size_t i = 0; i < deconv_ms_ptr->size(); i++){
    DeconvPeakPtr deconv_peak_ptr = deconv_ms_ptr->getPeakPtr(i);
    if(deconv_peak_ptr->getMonoMass() <= ext_min_mass) {
      double orig_mass = deconv_peak_ptr->getMonoMass();
      for (int j = 1; j < k; j++) {
        if (deconv_peak_ptr->getMonoMass() > new_prec_mass * j / k) {
          orig_mass -= sp_para_ptr->mod_mass_[j - 1];
        } 
      }
      ExtendPeakPtr extend_peak_ptr 
          = std::make_shared<ExtendPeak>(deconv_peak_ptr, orig_mass, 1.0);
      list.push_back(extend_peak_ptr);
    } else{
      for(size_t j = 0;j < ext_offsets.size();j++){
        double mass = deconv_peak_ptr->getMonoMass() + ext_offsets[j];
        for (int j = 1; j < k; j++) {
          if (deconv_peak_ptr->getMonoMass() + ext_offsets[j] > new_prec_mass * j / k) {
            mass -= sp_para_ptr->mod_mass_[j - 1];
          } 
        }
        ExtendPeakPtr extend_peak_ptr 
            = std::make_shared<ExtendPeak>(deconv_peak_ptr, mass, 1.0);
        list.push_back(extend_peak_ptr);
      }
    }
  }

  //filter extend_peak
  ExtendPeakPtrVec list_filtered;
  double min_mass = sp_para_ptr->getMinMass();
  double prec_mono_mass = header_ptr->getPrecMonoMass();
  for(size_t i =0; i < list.size();i++){
    double mass = list[i]->getPosition();
    if(mass >= min_mass && mass <= prec_mono_mass - min_mass){
      list_filtered.push_back(list[i]);
    }
  }

  // sort 
  std::sort(list_filtered.begin(),list_filtered.end(),ExtendPeak::cmpPosIncrease);

  //set error tolerance
  PeakTolerancePtr peak_tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  for (size_t i = 0; i < list_filtered.size();i++){
    double mass = list_filtered[i]->getBasePeakPtr()->getMonoMass();
    double ori_tole = peak_tole_ptr->compStrictErrorTole(mass);
    list_filtered[i]->setOrigTolerance(ori_tole);
    double reve_tole 
        = peak_tole_ptr->compRelaxErrorTole(mass, prec_mono_mass);
    list_filtered[i]->setReverseTolerance(reve_tole);
  }
  return ExtendMsPtr(new Ms<ExtendPeakPtr>(header_ptr,list_filtered));
}

ExtendMsPtrVec ExtendMsFactory::geneMsThreePtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                                  SpParaPtr sp_para_ptr, double new_prec_mass) {
  ExtendMsPtrVec extend_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    extend_ms_ptr_vec.push_back(geneMsThreePtr(deconv_ms_ptr_vec[i], sp_para_ptr, new_prec_mass));
  }
  return extend_ms_ptr_vec;
}

} /* namespace prot */
