// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#include <limits>

#include "base/logger.hpp"
#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "feature/feature.hpp"
#include "feature/feature_detect_mng.hpp"
#include "feature/feature_detect.hpp"
#include "feature/msalign_writer.hpp"

namespace prot {

void readSpectra(std::string &file_name, DeconvMsPtrVec &ms_ptr_vec) {
  int sp_num_in_group = 1;
  MsAlignReader sp_reader(file_name, sp_num_in_group, nullptr);

  DeconvMsPtr ms_ptr;
  //LOG_DEBUG("Start search");
  while((ms_ptr = sp_reader.getNextMs())!= nullptr){
    ms_ptr_vec.push_back(ms_ptr);
    //std::cout << std::flush <<  "reading spectrum " << header_ptr_vec.size() << "\r";
  }
  sp_reader.close();
  std::cout << std::endl;
}

void readHeaders(std::string &file_name, MsHeaderPtrVec &header_ptr_vec) {
  int sp_num_in_group = 1;
  MsAlignReader sp_reader(file_name, sp_num_in_group, nullptr);

  DeconvMsPtr ms_ptr;
  //LOG_DEBUG("Start search");
  while((ms_ptr = sp_reader.getNextMs())!= nullptr){
    header_ptr_vec.push_back(ms_ptr->getMsHeaderPtr());
    //std::cout << std::flush <<  "reading spectrum " << header_ptr_vec.size() << "\r";
  }
  sp_reader.close();
  std::cout << std::endl;
}

void outputHeaders(MsHeaderPtrVec &header_ptr_vec) {
  for (size_t i = 0; i < header_ptr_vec.size(); i++) {
    MsHeaderPtr ptr = header_ptr_vec[i];
    std::cout << ptr->getId() << "\t" << ptr->getFirstScanNum() << "\t";
    std::cout <<  ptr->getRetentionTime() << "\t" << ptr->getPrecMonoMass() << "\t";
    std::cout << ptr->getPrecMonoMz() << "\t" << ptr->getPrecCharge() << "\t" << ptr->getPrecInte() << std::endl; 
  }

}

bool isConsistent(MsHeaderPtr &a, MsHeaderPtr &b, FeatureDetectMngPtr mng_ptr) {
  double min_diff = std::numeric_limits<double>::max();
  std::vector<double> ext_masses = mng_ptr->getExtMasses(a->getPrecMonoMass());
  for (size_t i = 0; i < ext_masses.size(); i++) {
    double mass_diff = std::abs(ext_masses[i] - b->getPrecMonoMass());
    if (mass_diff < min_diff) {
      min_diff = mass_diff;
    }
  }

  double error_tole = mng_ptr->peak_tolerance_ptr_->compStrictErrorTole(a->getPrecMonoMass());
  if (min_diff <= error_tole) {
    return true;
  }
  return false;
}

bool containPrecursor(DeconvMsPtr ms1_ptr, MsHeaderPtr best_ptr, FeatureDetectMngPtr mng_ptr) {
  if (ms1_ptr == nullptr) return false;
  double prec_mass = best_ptr->getPrecMonoMass();
  //double prec_chrg = best_ptr->getPrecCharge();
  std::vector<double> ext_masses = mng_ptr->getExtMasses(prec_mass);

  double min_diff = std::numeric_limits<double>::max();
  for (size_t i = 0; i < ms1_ptr->size(); i++) {
    DeconvPeakPtr peak = ms1_ptr->getPeakPtr(i);
    //if (peak->getCharge() == prec_chrg) {
    // do not test charge 
    for (size_t j = 0; j < ext_masses.size(); j++) {
      double mass_diff = std::abs(ext_masses[j] - peak->getPosition());
      if (mass_diff < min_diff) {
        min_diff = mass_diff;
      }
    }
    //}
  }

  double error_tole = mng_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  if (min_diff <= error_tole) {
    return true;
  }
  return false;
}

int getMs1IdBegin(DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtr best_ptr, 
                  FeatureDetectMngPtr mng_ptr) {
  int cur_id = best_ptr->getMsOneId();
  while (cur_id > 0) {
    if (containPrecursor(ms1_ptr_vec[cur_id-1], best_ptr, mng_ptr)) {
      cur_id--;
    }
    else {
      break;
    }
  }
  return cur_id;
}

int getMs1IdEnd(DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtr best_ptr, 
                FeatureDetectMngPtr mng_ptr) {
  int cur_id = best_ptr->getMsOneId();
  while (cur_id < (int)ms1_ptr_vec.size()-1) {
    if (containPrecursor(ms1_ptr_vec[cur_id+1], best_ptr, mng_ptr)) {
      cur_id++;
    }
    else {
      break;
    }
  }
  return cur_id;
}

double getFeatureIntensity(DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtr best_ptr, 
                           int ms1_id_begin, int ms1_id_end, FeatureDetectMngPtr mng_ptr) {
  double sum = 0;
  double prec_mass = best_ptr->getPrecMonoMass();
  double error_tole = mng_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  std::vector<double> ext_masses = mng_ptr->getExtMasses(prec_mass);
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

int getMs2IdBegin(MsHeaderPtrVec &header_ptr_vec, MsHeaderPtr best_ptr, 
                  int ms1_id_begin) {
  int cur_id = best_ptr->getId();
  while (cur_id > 0) {
    if (header_ptr_vec[cur_id - 1]->getMsOneId() >= ms1_id_begin) {
      cur_id--;
    }
    else {
      break;
    }
  }
  return cur_id;
}

int getMs2IdEnd(MsHeaderPtrVec &header_ptr_vec, MsHeaderPtr best_ptr, int ms1_id_end) {
  int cur_id = best_ptr->getId();
  while (cur_id < (int)header_ptr_vec.size()-1) {
    if (header_ptr_vec[cur_id + 1]->getMsOneId() <= ms1_id_end) {
      cur_id++;
    }
    else {
      break;
    }
  }
  return cur_id;
}

void groupHeaders(DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtrVec &header_ptr_vec, 
                  FeatureDetectMngPtr mng_ptr, MsHeaderPtr2D &result_groups,
                  FeaturePtrVec &features) {
  MsHeaderPtrVec remain_ptrs = header_ptr_vec; 
  MsHeaderPtrVec sorted_ptrs = remain_ptrs; 
  while (sorted_ptrs.size() > 0) {
    //LOG_DEBUG("grouping");
    std::sort(sorted_ptrs.begin(), sorted_ptrs.end(), MsHeader::cmpPrecInteDec);
    MsHeaderPtr best_ptr = sorted_ptrs[0];
    MsHeaderPtrVec cur_group;
    cur_group.push_back(best_ptr);
    int best_id = best_ptr->getId();
    remain_ptrs[best_id] = nullptr;
    int ms1_id_begin = getMs1IdBegin(ms1_ptr_vec, best_ptr, mng_ptr);
    int ms1_id_end = getMs1IdEnd(ms1_ptr_vec, best_ptr, mng_ptr);
    double cur_inte = getFeatureIntensity(ms1_ptr_vec, best_ptr, ms1_id_begin,
                                          ms1_id_end, mng_ptr);

    int ms2_id_begin = getMs2IdBegin(header_ptr_vec, best_ptr, ms1_id_begin);
    int ms2_id_end = getMs2IdEnd(header_ptr_vec, best_ptr, ms1_id_end);

    for (int i = ms2_id_begin; i <= ms2_id_end; i++) {
      if (remain_ptrs[i] != nullptr) {
        if (isConsistent(best_ptr, remain_ptrs[i], mng_ptr)) {
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

    FeaturePtr feature_ptr(new Feature(feature_id, cur_mono_mass, cur_inte, ms1_scan_begin, ms1_scan_end));
    features.push_back(feature_ptr);
    sorted_ptrs.clear();
    for (size_t i = 0; i < remain_ptrs.size(); i++) {
      if (remain_ptrs[i] != nullptr) {
        sorted_ptrs.push_back(remain_ptrs[i]);
      }
    }
    //LOG_DEBUG("sorted ptr size " << sorted_ptrs.size());
  }
  /*for (size_t i = 0; i < result_groups.size(); i++) {*/
  //std::cout << "Group " << i << " number " << result_groups[i].size() << " ms1 scan begin " << features[i]->getScanBegin();
  //std::cout << " ms1 scan end " << features[i]->getScanEnd() << " inte " << features[i]->getIntensity() << std::endl;
  //for (size_t j = 0; j < result_groups[i].size(); j++) {
  //MsHeaderPtr ptr = result_groups[i][j];
  //std::cout << "\t" << ptr->getId() << "\t" << ptr->getFirstScanNum() << "\t";
  //std::cout <<  (ptr->getRetentionTime()/60) << "\t" << ptr->getPrecMonoMass() << "\t";
  //std::cout << ptr->getPrecMonoMz() << "\t" << ptr->getPrecCharge() << "\t" << ptr->getPrecInte() << std::endl; 
  //}
  /*}*/
}

void setFeatures(MsHeaderPtr2D &header_groups, FeaturePtrVec &features) {
  for (size_t i = 0; i < header_groups.size(); i++) {
    for (size_t j = 0; j < header_groups[i].size(); j++) {
      header_groups[i][j]->setFeatureId(features[i]->getId());
      header_groups[i][j]->setFeatureInte(features[i]->getIntensity());
    }
  }
}

void writeMsalign(std::string &input_file_name, std::string &output_file_name, MsHeaderPtrVec &header_ptrs) {
  int sp_num_in_group = 1;
  MsAlignReader sp_reader(input_file_name, sp_num_in_group, nullptr);
  std::ofstream of(output_file_name, std::ofstream::out);
  of.precision(16);
  DeconvMsPtr ms_ptr;
  int cnt = 0;
  while((ms_ptr = sp_reader.getNextMs())!= nullptr){
    ms_ptr->setHeaderPtr(header_ptrs[cnt]);
    MsalignWriter::writeText(of, ms_ptr);    
    cnt++;
  }
  sp_reader.close();
  of.close();
}

void FeatureDetect::process(std::string &xml_file_name){
  FeatureDetectMngPtr mng_ptr(new FeatureDetectMng());
  std::string base_name = FileUtil::basename(xml_file_name);
  // read ms1 deconvoluted spectra
  std::string ms1_file_name = base_name + ".ms1";
  DeconvMsPtrVec ms1_ptr_vec;
  readSpectra(ms1_file_name, ms1_ptr_vec);
  // read ms2 deconvoluted header
  LOG_DEBUG("start reading ms2");
  std::string ms2_file_name = base_name + ".ms2";
  MsHeaderPtrVec header_ptr_vec;
  readHeaders(ms2_file_name, header_ptr_vec);
  LOG_DEBUG("start grouping");
  MsHeaderPtr2D header_groups;
  FeaturePtrVec features;
  groupHeaders(ms1_ptr_vec, header_ptr_vec, mng_ptr, header_groups, features);
  setFeatures(header_groups, features);
  std::string output_file_name = base_name + ".msalign";
  writeMsalign(ms2_file_name, output_file_name, header_ptr_vec);
  //outputHeaders(header_ptr_vec);
}

//std::string output_file_name = FileUtil::basename(sp_file_name)+".feature";
}
