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

#include <limits>
#include <cmath>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/time_util.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_writer.hpp"
#include "deconv/feature/feature.hpp"
#include "deconv/feature/feature_para.hpp"
#include "deconv/feature/feature_detect.hpp"

namespace toppic {

namespace feature_detect {

void readSpectra(const std::string & file_name, DeconvMsPtrVec &ms_ptr_vec) {
  int sp_num_in_group = 1;
  MsAlignReader sp_reader(file_name, sp_num_in_group,
                          nullptr, std::set<std::string>());

  DeconvMsPtr ms_ptr;
  // LOG_DEBUG("Start search");
  while ((ms_ptr = sp_reader.getNextMs())!= nullptr) {
    ms_ptr->getMsHeaderPtr()->setMsLevel(1);
    ms_ptr_vec.push_back(ms_ptr);
    // std::cout << std::flush <<  "reading spectrum " << header_ptr_vec.size() << "\r";
  }
  sp_reader.close();
  // std::cout << std::endl;
}

void readHeaders(const std::string & file_name, MsHeaderPtrVec &header_ptr_vec) {
  int sp_num_in_group = 1;
  MsAlignReader sp_reader(file_name, sp_num_in_group, nullptr,
                          std::set<std::string>());

  DeconvMsPtr ms_ptr;
  LOG_DEBUG("Start search");
  while ((ms_ptr = sp_reader.getNextMs()) != nullptr) {
    header_ptr_vec.push_back(ms_ptr->getMsHeaderPtr());
    std::cout << std::flush <<  "reading spectrum " << header_ptr_vec.size() << "\r";
  }
  sp_reader.close();
  std::cout << std::endl;
}

void outputHeaders(const MsHeaderPtrVec &header_ptr_vec) {
  for (size_t i = 0; i < header_ptr_vec.size(); i++) {
    MsHeaderPtr ptr = header_ptr_vec[i];
    std::cout << ptr->getId() << "\t" << ptr->getFirstScanNum() << "\t";
    std::cout <<  ptr->getRetentionTime() << "\t" << ptr->getPrecMonoMass() << "\t";
    std::cout << ptr->getPrecMonoMz() << "\t" << ptr->getPrecCharge() << "\t"
        << ptr->getPrecInte() << std::endl;
  }
}

bool isConsistent(MsHeaderPtr a, MsHeaderPtr b, FeatureParaPtr para_ptr) {
  double min_diff = std::numeric_limits<double>::max();
  std::vector<double> ext_masses = para_ptr->getExtMasses(a->getPrecMonoMass());
  for (size_t i = 0; i < ext_masses.size(); i++) {
    double mass_diff = std::abs(ext_masses[i] - b->getPrecMonoMass());
    if (mass_diff < min_diff) {
      min_diff = mass_diff;
    }
  }

  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(a->getPrecMonoMass());
  if (min_diff <= error_tole) {
    return true;
  }
  return false;
}

bool containPrecursor(DeconvMsPtr ms1_ptr, MsHeaderPtr best_ptr, FeatureParaPtr para_ptr) {
  if (ms1_ptr == nullptr) return false;
  double prec_mass = best_ptr->getPrecMonoMass();
  // double prec_chrg = best_ptr->getPrecCharge();
  std::vector<double> ext_masses = para_ptr->getExtMasses(prec_mass);

  double min_diff = std::numeric_limits<double>::max();
  for (size_t i = 0; i < ms1_ptr->size(); i++) {
    DeconvPeakPtr peak = ms1_ptr->getPeakPtr(i);
    // if (peak->getCharge() == prec_chrg) {
    // do not test charge
    for (size_t j = 0; j < ext_masses.size(); j++) {
      double mass_diff = std::abs(ext_masses[j] - peak->getPosition());
      if (mass_diff < min_diff) {
        min_diff = mass_diff;
      }
    }
    //}
  }

  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  if (min_diff <= error_tole) {
    return true;
  }
  return false;
}

int getMs1IdBegin(const DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtr best_ptr,
                  FeatureParaPtr para_ptr) {
  int cur_id = best_ptr->getMsOneId();
  while (cur_id > 0) {
    if (containPrecursor(ms1_ptr_vec[cur_id-1], best_ptr, para_ptr)) {
      cur_id--;
    } else {
      break;
    }
  }
  return cur_id;
}

int getMs1IdEnd(const DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtr best_ptr,
                FeatureParaPtr para_ptr) {
  int cur_id = best_ptr->getMsOneId();
  while (cur_id < static_cast<int>(ms1_ptr_vec.size()) - 1) {
    if (containPrecursor(ms1_ptr_vec[cur_id + 1], best_ptr, para_ptr)) {
      cur_id++;
    } else {
      break;
    }
  }
  return cur_id;
}

double getFeatureIntensity(const DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtr best_ptr,
                           int ms1_id_begin, int ms1_id_end, FeatureParaPtr para_ptr) {
  if (ms1_ptr_vec.size() == 0) return 0.0;
  double sum = 0;
  double prec_mass = best_ptr->getPrecMonoMass();
  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  std::vector<double> ext_masses = para_ptr->getExtMasses(prec_mass);
  for (int i = ms1_id_begin; i <= ms1_id_end; i++) {
    DeconvMsPtr ms1_ptr = ms1_ptr_vec[i];
    for (size_t j = 0; j < ms1_ptr->size(); j++) {
      DeconvPeakPtr peak = ms1_ptr->getPeakPtr(j);
      for (size_t k = 0; k < ext_masses.size(); k++) {
        double mass_diff = std::abs(ext_masses[k] - peak->getPosition());
        if (mass_diff <= error_tole) {
          sum += peak->getIntensity();
          break;
        }
      }
    }
  }
  return sum;
}

int getMs2IdBegin(const MsHeaderPtrVec &header_ptr_vec, MsHeaderPtr best_ptr,
                  int ms1_id_begin) {
  int cur_id = best_ptr->getId();
  while (cur_id > 0) {
    if (header_ptr_vec[cur_id - 1]->getMsOneId() >= ms1_id_begin) {
      cur_id--;
    } else {
      break;
    }
  }
  return cur_id;
}

int getMs2IdEnd(MsHeaderPtrVec &header_ptr_vec, MsHeaderPtr best_ptr, int ms1_id_end) {
  int cur_id = best_ptr->getId();
  while (cur_id < static_cast<int>(header_ptr_vec.size()) - 1) {
    if (header_ptr_vec[cur_id + 1]->getMsOneId() <= ms1_id_end) {
      cur_id++;
    } else {
      break;
    }
  }
  return cur_id;
}

void groupHeaders(DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtrVec &header_ptr_vec,
                  FeatureParaPtr para_ptr, MsHeaderPtr2D &result_groups,
                  FeaturePtrVec &features) {
  if (ms1_ptr_vec.size() == 0) return;
  MsHeaderPtrVec remain_ptrs = header_ptr_vec;
  MsHeaderPtrVec sorted_ptrs = remain_ptrs;
  while (sorted_ptrs.size() > 0) {
    // LOG_DEBUG("grouping");
    std::sort(sorted_ptrs.begin(), sorted_ptrs.end(), MsHeader::cmpPrecInteDec);
    MsHeaderPtr best_ptr = sorted_ptrs[0];
    MsHeaderPtrVec cur_group;
    cur_group.push_back(best_ptr);
    int best_id = best_ptr->getId();
    remain_ptrs[best_id] = nullptr;
    int ms1_id_begin = getMs1IdBegin(ms1_ptr_vec, best_ptr, para_ptr);
    int ms1_id_end = getMs1IdEnd(ms1_ptr_vec, best_ptr, para_ptr);
    double cur_inte = getFeatureIntensity(ms1_ptr_vec, best_ptr, ms1_id_begin,
                                          ms1_id_end, para_ptr);
    if (cur_inte == 0) {
      int spec_id = ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr()->getId();
      DeconvPeakPtrVec peak_vec = ms1_ptr_vec[ms1_id_begin]->getPeakPtrVec();
      peak_vec.push_back(std::make_shared<DeconvPeak>(spec_id, 
                                                      peak_vec.size(),
                                                      best_ptr->getPrecMonoMass(),
                                                      best_ptr->getPrecInte(),
                                                      best_ptr->getPrecCharge()));
      ms1_ptr_vec[ms1_id_begin]
          = std::make_shared<Ms<DeconvPeakPtr> >(ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr(),
                                                 peak_vec);
      cur_inte = best_ptr->getPrecInte();
    }
    int ms2_id_begin = getMs2IdBegin(header_ptr_vec, best_ptr, ms1_id_begin);
    int ms2_id_end = getMs2IdEnd(header_ptr_vec, best_ptr, ms1_id_end);

    for (int i = ms2_id_begin; i <= ms2_id_end; i++) {
      if (remain_ptrs[i] != nullptr) {
        if (isConsistent(best_ptr, remain_ptrs[i], para_ptr)) {
          cur_group.push_back(remain_ptrs[i]);
          remain_ptrs[i] = nullptr;
        }
      }
    }

    result_groups.push_back(cur_group);
    int feature_id = features.size();
    double cur_mono_mass = best_ptr->getPrecMonoMass();
    int ms1_scan_begin = ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr()->getFirstScanNum();
    int ms1_scan_end = ms1_ptr_vec[ms1_id_end]->getMsHeaderPtr()->getFirstScanNum();

    FeaturePtr feature_ptr = std::make_shared<Feature>(feature_id, cur_mono_mass, cur_inte,
                                                       ms1_scan_begin, ms1_scan_end);
    features.push_back(feature_ptr);
    sorted_ptrs.clear();
    for (size_t i = 0; i < remain_ptrs.size(); i++) {
      if (remain_ptrs[i] != nullptr) {
        sorted_ptrs.push_back(remain_ptrs[i]);
      }
    }
    // LOG_DEBUG("sorted ptr size " << sorted_ptrs.size());
  }
  /*for (size_t i = 0; i < result_groups.size(); i++) {
  //std::cout << "Group " << i << " number " << result_groups[i].size() << " ms1 scan begin " << features[i]->getScanBegin();
  //std::cout << " ms1 scan end " << features[i]->getScanEnd() << " inte " << features[i]->getIntensity() << std::endl;
  //for (size_t j = 0; j < result_groups[i].size(); j++) {
  //MsHeaderPtr ptr = result_groups[i][j];
  //std::cout << "\t" << ptr->getId() << "\t" << ptr->getFirstScanNum() << "\t";
  //std::cout <<  (ptr->getRetentionTime()/60) << "\t" << ptr->getPrecMonoMass() << "\t";
  //std::cout << ptr->getPrecMonoMz() << "\t" << ptr->getPrecCharge() << "\t" << ptr->getPrecInte() << std::endl; 
  //}
  }*/
}

void setFeatures(MsHeaderPtr2D &header_groups, const FeaturePtrVec &features) {
  for (size_t i = 0; i < header_groups.size(); i++) {
    for (size_t j = 0; j < header_groups[i].size(); j++) {
      header_groups[i][j]->setFeatureId(features[i]->getId());
      header_groups[i][j]->setFeatureInte(features[i]->getIntensity());
    }
  }
}

void writeMSFT(const std::string & input_file_name,
               const std::string & output_file_name,
               const MsHeaderPtrVec &header_ptrs,
               std::string argu_str) {
  int sp_num_in_group = 1;
  MsAlignReader sp_reader(input_file_name, sp_num_in_group, nullptr, std::set<std::string>());
  std::ofstream of(output_file_name, std::ofstream::out);
  of.precision(16);
  time_util::addTimeStamp(argu_str);
  of << argu_str;
  of << "ID" << "\t"
      << "SCANS" << "\t"
      << "MS_ONE_ID" << "\t"
      << "MS_ONE_SCAN" << "\t"
      << "PRECURSOR_MASS" << "\t"
      << "PRECURSOR_INTENSITY" << "\t"
      << "FEATURE_ID" << "\t"
      << "FEATURE_INTENSITY" << std::endl;
  DeconvMsPtr ms_ptr;
  int cnt = 0;
  while ((ms_ptr = sp_reader.getNextMs()) != nullptr) {
    ms_ptr->setHeaderPtr(header_ptrs[cnt]);
    of << header_ptrs[cnt]->getId() << "\t"
        << header_ptrs[cnt]->getScansString() << "\t";
    if (header_ptrs[cnt]->getMsLevel() > 1) {
      of << header_ptrs[cnt]->getMsOneId() << "\t"
          << header_ptrs[cnt]->getMsOneScan() << "\t"
          <<  header_ptrs[cnt]->getPrecMonoMass() << "\t"
          << header_ptrs[cnt]->getPrecInte() << "\t";
      if (header_ptrs[cnt]->getFeatureId() >= 0) {
        of << header_ptrs[cnt]->getFeatureId() << "\t"
            << header_ptrs[cnt]->getFeatureInte() << std::endl;
      } else {
        of << "-" << "\t" << "-" << std::endl;
      }
    } else {
      of << "-" << "\t"
          << "-" << "\t"
          << "-" << "\t"
          << "-" << "\t"
          << "-" << "\t"
          << "-" << std::endl;
    }
    cnt++;
  }
  sp_reader.close();
  of.close();
}

void process(std::string &sp_file_name, bool missing_level_one, 
             std::string &argu_str) {
  //logger::setLogLevel(2);
  FeatureParaPtr para_ptr = std::make_shared<FeaturePara>();
  std::string base_name = file_util::basename(sp_file_name);
  // read ms1 deconvoluted spectra
  std::string ms1_file_name = base_name + "_ms1.msalign";
  DeconvMsPtrVec ms1_ptr_vec;
  if (!missing_level_one) readSpectra(ms1_file_name, ms1_ptr_vec);
  // read ms2 deconvoluted header
  LOG_DEBUG("start reading ms2");
  std::string ms2_file_name = base_name + "_ms2.msalign";
  MsHeaderPtrVec header_ptr_vec;
  readHeaders(ms2_file_name, header_ptr_vec);
  LOG_DEBUG("start grouping");
  MsHeaderPtr2D header_groups;
  FeaturePtrVec features;
  groupHeaders(ms1_ptr_vec, header_ptr_vec, para_ptr, header_groups, features);
  setFeatures(header_groups, features);
  std::string output_file_name = base_name + ".feature";
  writeMSFT(ms2_file_name, output_file_name, header_ptr_vec, argu_str);
  std::ofstream of1(ms1_file_name, std::ofstream::out);
  of1.precision(16);
  time_util::addTimeStamp(argu_str);
  of1 << argu_str;
  for (size_t i = 0; i < ms1_ptr_vec.size(); i++) {
    msalign_writer::write(of1, ms1_ptr_vec[i]);
  }
  of1.close();
  // outputHeaders(header_ptr_vec);
}

}

}  // namespace toppic
