//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <numeric>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/base/mass_constant.hpp"
#include "ms/spec/base_peak_type.hpp"
#include "ms/factory/prm_ms_factory.hpp"

namespace toppic {

namespace prm_ms_factory {

void addTwoMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr,
                  double prec_mono_mass, double n_term_label_mass, 
                  ActivationPtr active_type_ptr, PeakTolerancePtr tole_ptr) {

  double peak_mass = deconv_peak_ptr->getMonoMass();
  double n_term_mass = peak_mass - active_type_ptr->getN_BYShift() - n_term_label_mass;
  double default_score = 1.0;
  PrmPeakPtr new_peak_ptr
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                  BasePeakType::ORIGINAL, n_term_mass, default_score);
  new_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(peak_mass));
  new_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(peak_mass));
  new_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(peak_mass, prec_mono_mass));
  list.push_back(new_peak_ptr);
  double reverse_mass = prec_mono_mass - (peak_mass -active_type_ptr->getC_BYShift()) - n_term_label_mass;
  PrmPeakPtr reverse_peak_ptr
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                  BasePeakType::REVERSED, reverse_mass, default_score);
  reverse_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(peak_mass));
  reverse_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(peak_mass));
  reverse_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(peak_mass, prec_mono_mass));
  list.push_back(reverse_peak_ptr);
}

void addSuffixTwoMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr,
                        double prec_mono_mass, ActivationPtr active_type_ptr,
                        PeakTolerancePtr tole_ptr) {
  double default_score = 1.0;
  double peak_mass = deconv_peak_ptr->getMonoMass();
  double c_res_mass = peak_mass - active_type_ptr->getC_BYShift() - mass_constant::getWaterMass();
  PrmPeakPtr new_peak_ptr
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr, BasePeakType::ORIGINAL, c_res_mass, default_score);
  new_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(peak_mass));
  new_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(peak_mass, prec_mono_mass));
  new_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(peak_mass));
  list.push_back(new_peak_ptr);
  double reverse_mass = prec_mono_mass - (peak_mass - active_type_ptr->getN_BYShift()) 
                        - mass_constant::getWaterMass();
  PrmPeakPtr reverse_peak_ptr
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                  BasePeakType::REVERSED, reverse_mass, default_score);
  reverse_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(peak_mass));
  reverse_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(peak_mass));
  reverse_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(peak_mass, prec_mono_mass));
  list.push_back(reverse_peak_ptr);
}

void addSixMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr,
                  double prec_mono_mass, double n_term_label_mass, 
                  ActivationPtr active_type_ptr, PeakTolerancePtr tole_ptr, 
                  const std::vector<double> &offsets) {
  double default_score = 1.0;
  double peak_mass = deconv_peak_ptr->getMonoMass();
  for (size_t i = 0; i < offsets.size(); i++) {
    double mass = peak_mass - active_type_ptr->getN_BYShift() - n_term_label_mass + offsets[i];
    PrmPeakPtr peak_ptr = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                                    BasePeakType::ORIGINAL, mass, default_score);
    peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(peak_mass));
    peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(peak_mass));
    peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(peak_mass, prec_mono_mass));
    list.push_back(peak_ptr);
  }
  for (size_t i = 0; i < offsets.size(); i++) {
    double mass = prec_mono_mass-(peak_mass-active_type_ptr->getC_BYShift() - n_term_label_mass +offsets[i]);
    PrmPeakPtr peak_ptr = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                                    BasePeakType::REVERSED, mass, default_score);
    peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(peak_mass));
    peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(peak_mass));
    peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(peak_mass, prec_mono_mass));
    list.push_back(peak_ptr);
  }
}

void addSuffixSixMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr,
                        double prec_mono_mass, ActivationPtr active_type_ptr,
                        PeakTolerancePtr tole_ptr, const std::vector<double> &offsets) {
  double default_score = 1.0;
  double ori_mass = deconv_peak_ptr->getMonoMass();
  double c_res_mass = ori_mass - (active_type_ptr->getC_BYShift() + mass_constant::getWaterMass());
  for (size_t i = 0; i < offsets.size(); i++) {
    double mass = c_res_mass + offsets[i];
    PrmPeakPtr peak_ptr = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                                    BasePeakType::ORIGINAL, mass, default_score);
    peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
    peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    list.push_back(peak_ptr);
  }
  double reverse_mass = prec_mono_mass - (ori_mass-active_type_ptr->getN_BYShift()) 
    - mass_constant::getWaterMass();
  for (size_t i = 0; i < offsets.size(); i++) {
    double mass = reverse_mass + offsets[i];
    PrmPeakPtr peak_ptr = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                                    BasePeakType::REVERSED, mass, default_score);
    peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
    peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    list.push_back(peak_ptr);
  }
}

void filterPeaks(const PrmPeakPtrVec &peak_list, PrmPeakPtrVec &filtered_list,
                 double prec_mono_mass, double min_mass) {
  for (size_t i = 0; i < peak_list.size(); i++) {
    if (peak_list[i]->getPosition() >= min_mass &&
        peak_list[i]->getPosition() <= prec_mono_mass - min_mass) {
      filtered_list.push_back(peak_list[i]);
    }
  }
}

void geneApproxSpecPeaks(PrmPeakPtrVec &peak_list, double prec_mono_mass, 
                         const std::vector<double> &mod_mass) {
  for (size_t i = 0; i < peak_list.size(); i++) {
    double mass = peak_list[i]->getMonoMass();
    if (peak_list[i]->getMonoMass() > prec_mono_mass * 1 / 6) {
      mass -= mod_mass[0];
    }
    if (peak_list[i]->getMonoMass() > prec_mono_mass * 3 / 6) {
      mass -= mod_mass[1];
    }
    if (peak_list[i]->getMonoMass() > prec_mono_mass * 5 / 6) {
      mass -= mod_mass[2];
    }
    peak_list[i] = std::make_shared<PrmPeak>(peak_list[i]->getSpectrumId(),
                                             peak_list[i]->getBasePeakPtr(),
                                             peak_list[i]->getBaseTypePtr(),
                                             mass,
                                             peak_list[i]->getScore(),
                                             peak_list[i]->getStrictTolerance(),
                                             peak_list[i]->getNStrictCRelaxTolerance(),
                                             peak_list[i]->getNRelaxCStrictTolerance());
  }
}

PrmMsPtr geneMsTwoPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                      double prec_mono_mass, double n_term_label_mass, 
                      const std::vector<double> & mod_mass) {
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  double prec_mass_without_label = prec_mono_mass - n_term_label_mass;
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, prec_mass_without_label);
  // get two prm peaks 
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  PrmPeakPtrVec list;
  for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
    addTwoMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                 n_term_label_mass, active_type_ptr, tole_ptr);
  }
  // filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mass_without_label, sp_para_ptr->getMinMass());
  std::sort(list_filtered.begin(), list_filtered.end(), PrmPeak::cmpPosInc);
  if (mod_mass.size() > 0) {
    // generate approximate spectrum
    geneApproxSpecPeaks(list_filtered, prec_mass_without_label, mod_mass);
  }
  return std::make_shared<Ms<PrmPeakPtr> >(header_ptr, list_filtered);
}

PrmMsPtr geneSuffixMsTwoPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                            double prec_mono_mass, double n_term_label_mass, 
                            const std::vector<double> & mod_mass) {
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  double prec_mass_without_label = prec_mono_mass - n_term_label_mass;
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, prec_mass_without_label);
  // get two prm peaks 
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  PrmPeakPtrVec list;
  for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
    addSuffixTwoMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                       active_type_ptr, tole_ptr);
  }
  // filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mass_without_label, sp_para_ptr->getMinMass());
  std::sort(list_filtered.begin(), list_filtered.end(), PrmPeak::cmpPosInc);
  if (mod_mass.size() > 0) {
    // generate approximate spectrum using mod masses in the reversed order
    std::vector<double> rev_mod_mass = mod_mass;
    std::reverse(rev_mod_mass.begin(), rev_mod_mass.end()); 
    geneApproxSpecPeaks(list_filtered, prec_mass_without_label, rev_mod_mass);
  }
  return std::make_shared<Ms<PrmPeakPtr> >(header_ptr, list_filtered);
}

PrmMsPtr geneMsSixPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                      double prec_mono_mass, double n_term_label_mass) {
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  double prec_mass_without_label = prec_mono_mass - n_term_label_mass;
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, prec_mass_without_label);
  // get six prm peaks 
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  double extend_min_mass = sp_para_ptr->getExtendMinMass();
  PrmPeakPtrVec list;
  for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
    if (deconv_ms_ptr->getPeakPtr(i)->getMonoMass() <= extend_min_mass) {
      addTwoMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                   n_term_label_mass, active_type_ptr, tole_ptr);
    } else {
      addSixMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                   n_term_label_mass, active_type_ptr, tole_ptr, sp_para_ptr->getExtendOffsets());
    }
  }

  // filter peaks 
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mass_without_label, sp_para_ptr->getMinMass());
  std::sort(list_filtered.begin(), list_filtered.end(), PrmPeak::cmpPosInc);

  return std::make_shared<Ms<PrmPeakPtr>>(header_ptr, list_filtered);
}

PrmMsPtr geneSuffixMsSixPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                            double prec_mono_mass, double n_term_label_mass) {
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  double prec_mass_without_label = prec_mono_mass - n_term_label_mass;
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, prec_mass_without_label);
  // get two prm peaks 
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  double extend_min_mass = sp_para_ptr->getExtendMinMass();
  PrmPeakPtrVec list;
  for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
    if (deconv_ms_ptr->getPeakPtr(i)->getMonoMass() <= extend_min_mass) {
      addSuffixTwoMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                         active_type_ptr, tole_ptr);
    } else {
      addSuffixSixMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                         active_type_ptr, tole_ptr, sp_para_ptr->getExtendOffsets());
    }
  }
  // filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mass_without_label, sp_para_ptr->getMinMass());
  std::sort(list_filtered.begin(), list_filtered.end(), PrmPeak::cmpPosInc);
  return std::make_shared<Ms<PrmPeakPtr> >(header_ptr, list_filtered);
}


PrmMsPtrVec geneMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                            SpParaPtr sp_para_ptr, double prec_mono_mass,
                            double n_term_label_mass, const std::vector<double> & mod_mass) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneMsTwoPtr(deconv_ms_ptr_vec[i], i, sp_para_ptr, 
                                          prec_mono_mass, n_term_label_mass, mod_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec geneSuffixMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                                  SpParaPtr sp_para_ptr, double prec_mono_mass,
                                  double n_term_label_mass,
                                  const std::vector<double> & mod_mass) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneSuffixMsTwoPtr(deconv_ms_ptr_vec[i], i,
                                                sp_para_ptr, prec_mono_mass, 
                                                n_term_label_mass, mod_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec geneMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                            SpParaPtr sp_para_ptr, double prec_mono_mass,
                            double n_term_label_mass) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneMsSixPtr(deconv_ms_ptr_vec[i], i, sp_para_ptr, 
                                          prec_mono_mass, n_term_label_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec geneSuffixMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                                  SpParaPtr sp_para_ptr, double prec_mono_mass, 
                                  double n_term_label_mass) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneSuffixMsSixPtr(deconv_ms_ptr_vec[i], i,
                                                sp_para_ptr, prec_mono_mass, n_term_label_mass));
  }
  return prm_ms_ptr_vec;
}

}  // namespace prm_ms_factory
}  // namespace toppic
