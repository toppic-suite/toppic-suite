//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include <algorithm>

#include "common/util/logger.hpp"
#include "common/base/mass_constant.hpp"
#include "spec/base_peak_type.hpp"
#include "spec/prm_ms_factory.hpp"

namespace toppic {

namespace prm_ms_factory {

void addTwoMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr,
                  double prec_mono_mass, ActivationPtr active_type_ptr, PeakTolerancePtr tole_ptr) {
  double ori_mass = deconv_peak_ptr->getMonoMass();
  double n_term_mass = ori_mass - active_type_ptr->getNShift();
  PrmPeakPtr new_peak_ptr
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                  BasePeakType::ORIGINAL, n_term_mass, 1);
  new_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  new_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  new_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  list.push_back(new_peak_ptr);
  double reverse_mass = prec_mono_mass - (deconv_peak_ptr->getMonoMass()-active_type_ptr->getCShift());
  PrmPeakPtr reverse_peak_ptr
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                  BasePeakType::REVERSED, reverse_mass, 1);
  reverse_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  list.push_back(reverse_peak_ptr);
}

void addSuffixTwoMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr,
                        double prec_mono_mass, ActivationPtr active_type_ptr,
                        PeakTolerancePtr tole_ptr) {
  double ori_mass = deconv_peak_ptr->getMonoMass();
  double c_res_mass = ori_mass - active_type_ptr->getCShift() - mass_constant::getWaterMass();
  PrmPeakPtr new_peak_ptr
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr, BasePeakType::ORIGINAL, c_res_mass, 1);
  new_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  new_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  new_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  list.push_back(new_peak_ptr);
  double reverse_mass = prec_mono_mass - (deconv_peak_ptr->getMonoMass()-active_type_ptr->getNShift())
      - mass_constant::getWaterMass();
  PrmPeakPtr reverse_peak_ptr
      = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                  BasePeakType::REVERSED, reverse_mass, 1);
  reverse_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(ori_mass));
  reverse_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
  list.push_back(reverse_peak_ptr);
}

void addSixMasses(PrmPeakPtrVec &list, int spec_id, DeconvPeakPtr deconv_peak_ptr,
                  double prec_mono_mass, ActivationPtr active_type_ptr,
                  PeakTolerancePtr tole_ptr, const std::vector<double> &offsets) {
  double ori_mass = deconv_peak_ptr->getMonoMass();
  for (size_t i = 0; i < offsets.size(); i++) {
    double mass = ori_mass - active_type_ptr->getNShift() + offsets[i];
    PrmPeakPtr peak_ptr = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                                    BasePeakType::ORIGINAL, mass, 1);
    peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
    list.push_back(peak_ptr);
  }
  for (size_t i = 0; i < offsets.size(); i++) {
    double mass = prec_mono_mass-(ori_mass-active_type_ptr->getCShift()+offsets[i]);
    PrmPeakPtr peak_ptr = std::make_shared<PrmPeak>(spec_id, deconv_peak_ptr,
                                                    BasePeakType::REVERSED, mass, 1);
    peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(ori_mass));
    peak_ptr->setNStrictCRelacTolerance(tole_ptr->compRelaxErrorTole(ori_mass, prec_mono_mass));
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

PrmMsPtr geneMsTwoPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                      double prec_mono_mass, const std::vector<double> & mod_mass) {
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  double new_prec_mono_mass = prec_mono_mass - std::accumulate(mod_mass.begin(), mod_mass.end(), 0.0);
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, new_prec_mono_mass);
  // getSpTwoPrmPeak
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  PrmPeakPtrVec list;
  for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
    addTwoMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                 active_type_ptr, tole_ptr);
  }
  // filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());
  std::sort(list_filtered.begin(), list_filtered.end(), PrmPeak::cmpPosInc);
  if (mod_mass.size() > 0) {
    for (size_t i = 0; i < list_filtered.size(); i++) {
      double mass = list_filtered[i]->getMonoMass();

      if (list_filtered[i]->getMonoMass() > prec_mono_mass / 1) {
        mass -= mod_mass[0];
      }

      if (list_filtered[i]->getMonoMass() > prec_mono_mass * 3 / 6) {
        mass -= mod_mass[1];
      }

      if (list_filtered[i]->getMonoMass() > prec_mono_mass * 5 / 6) {
        mass -= mod_mass[2];
      }

      list_filtered[i] = std::make_shared<PrmPeak>(list_filtered[i]->getSpectrumId(),
                                                   list_filtered[i]->getBasePeakPtr(),
                                                   list_filtered[i]->getBaseTypePtr(),
                                                   mass,
                                                   list_filtered[i]->getScore(),
                                                   list_filtered[i]->getStrictTolerance(),
                                                   list_filtered[i]->getNStrictCRelaxTolerance(),
                                                   list_filtered[i]->getNRelaxCStrictTolerance());
    }
  }
  return std::make_shared<Ms<PrmPeakPtr> >(header_ptr, list_filtered);
}

PrmMsPtr geneSuffixMsTwoPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                            double prec_mono_mass, const std::vector<double> & mod_mass) {
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  double new_prec_mono_mass = prec_mono_mass - std::accumulate(mod_mass.begin(), mod_mass.end(), 0.0);
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, new_prec_mono_mass);
  // getSpTwoPrmPeak
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  PrmPeakPtrVec list;
  for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
    addSuffixTwoMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), new_prec_mono_mass,
                       active_type_ptr, tole_ptr);
  }
  // filter low mass peaks
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());
  std::sort(list_filtered.begin(), list_filtered.end(), PrmPeak::cmpPosInc);
  if (mod_mass.size() > 0) {
    for (size_t i = 0; i < list_filtered.size(); i++) {
      double mass = list_filtered[i]->getMonoMass();

      if (list_filtered[i]->getMonoMass() > prec_mono_mass * 1 / 6) {
        mass -= mod_mass[2];
      }

      if (list_filtered[i]->getMonoMass() > prec_mono_mass * 3 / 6) {
        mass -= mod_mass[1];
      }

      if (list_filtered[i]->getMonoMass() > prec_mono_mass * 5 / 6) {
        mass -= mod_mass[0];
      }

      list_filtered[i] = std::make_shared<PrmPeak>(list_filtered[i]->getSpectrumId(),
                                                   list_filtered[i]->getBasePeakPtr(),
                                                   list_filtered[i]->getBaseTypePtr(),
                                                   mass,
                                                   list_filtered[i]->getScore(),
                                                   list_filtered[i]->getStrictTolerance(),
                                                   list_filtered[i]->getNStrictCRelaxTolerance(),
                                                   list_filtered[i]->getNRelaxCStrictTolerance());
    }
  }
  return std::make_shared<Ms<PrmPeakPtr> >(header_ptr, list_filtered);
}

PrmMsPtr geneMsSixPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                      double prec_mono_mass) {
  MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
  MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, prec_mono_mass);
  // getSpSixPrmPeak
  ActivationPtr active_type_ptr = header_ptr->getActivationPtr();
  PeakTolerancePtr tole_ptr = sp_para_ptr->getPeakTolerancePtr();
  double extend_min_mass = sp_para_ptr->getExtendMinMass();
  PrmPeakPtrVec list;
  for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
    if (deconv_ms_ptr->getPeakPtr(i)->getMonoMass() <= extend_min_mass) {
      addTwoMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                   active_type_ptr, tole_ptr);
    } else {
      addSixMasses(list, spec_id, deconv_ms_ptr->getPeakPtr(i), prec_mono_mass,
                   active_type_ptr, tole_ptr, sp_para_ptr->getExtendOffsets());
    }
  }

  // filterPrmPeak
  PrmPeakPtrVec list_filtered;
  filterPeaks(list, list_filtered, prec_mono_mass, sp_para_ptr->getMinMass());
  std::sort(list_filtered.begin(), list_filtered.end(), PrmPeak::cmpPosInc);

  return std::make_shared<Ms<PrmPeakPtr>>(header_ptr, list_filtered);
}

PrmMsPtr geneShiftMsSixPtr(DeconvMsPtr deconv_ms_ptr, int spec_id, SpParaPtr sp_para_ptr,
                           double prec_mono_mass, double shift) {
  PrmMsPtr prm_ms_ptr = geneMsSixPtr(deconv_ms_ptr, spec_id, sp_para_ptr, prec_mono_mass);
  MsHeaderPtr ori_header_ptr = prm_ms_ptr->getMsHeaderPtr();
  MsHeaderPtr header_ptr = std::make_shared<MsHeader>(*ori_header_ptr.get());
  double mono_mz = (header_ptr->getPrecMonoMass()+shift) /header_ptr->getPrecCharge();
  header_ptr->setPrecMonoMz(mono_mz);
  PrmPeakPtrVec prm_peak_list;
  for (size_t i = 0; i < prm_ms_ptr->size(); i++) {
    double pos = prm_ms_ptr->getPeakPtr(i)->getPosition()+shift;
    if (pos > 0) {
      prm_ms_ptr->getPeakPtr(i)->setPosition(pos);
      prm_peak_list.push_back(prm_ms_ptr->getPeakPtr(i));
    }
  }
  return std::make_shared<Ms<PrmPeakPtr>>(header_ptr, prm_peak_list);
}

PrmMsPtrVec geneMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                            SpParaPtr sp_para_ptr, double prec_mono_mass,
                            const std::vector<double> & mod_mass) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneMsTwoPtr(deconv_ms_ptr_vec[i], i,
                                          sp_para_ptr, prec_mono_mass, mod_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec geneSuffixMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                                  SpParaPtr sp_para_ptr, double prec_mono_mass,
                                  const std::vector<double> & mod_mass) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneSuffixMsTwoPtr(deconv_ms_ptr_vec[i], i,
                                                sp_para_ptr, prec_mono_mass, mod_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec geneMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                            SpParaPtr sp_para_ptr, double prec_mono_mass) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneMsSixPtr(deconv_ms_ptr_vec[i], i,
                                          sp_para_ptr, prec_mono_mass));
  }
  return prm_ms_ptr_vec;
}

PrmMsPtrVec geneShiftMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                                 SpParaPtr sp_para_ptr, double prec_mono_mass,
                                 double shift) {
  PrmMsPtrVec prm_ms_ptr_vec;
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    prm_ms_ptr_vec.push_back(geneShiftMsSixPtr(deconv_ms_ptr_vec[i], i, sp_para_ptr,
                                               prec_mono_mass, shift));
  }
  return prm_ms_ptr_vec;
}

}  // namespace prm_ms_factory
}  // namespace toppic
