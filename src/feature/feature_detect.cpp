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

#include <cmath>
#include <iostream>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/time_util.hpp"
#include "common/base/mod_util.hpp"
#include "seq/fasta_util.hpp"
#include "seq/fasta_index_reader.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/msalign_reader.hpp"
#include "feature/feature_para.hpp"
#include "feature/frac_feature.hpp"
#include "feature/frac_feature_writer.hpp"
#include "feature/spec_feature.hpp"
#include "feature/spec_feature_writer.hpp"
#include "feature/feature_detect.hpp"

namespace toppic {

namespace feature_detect {

bool peakExists (DeconvMsPtrVec &ms1_ptr_vec, DeconvPeakPtr peak) {
  int sp_id = peak->getSpId();
  int peak_id = peak->getId();
  DeconvPeakPtr remain_peak = ms1_ptr_vec[sp_id]->getPeakPtr(peak_id);
  if (remain_peak == nullptr) {
    return false;
  }
  else {
    return true;
  }
}

void getMatchedPeaks(DeconvMsPtrVec &ms1_ptr_vec, double prec_mass,
                     DeconvPeakPtrVec &matched_peaks, 
                     int ms1_id_begin, int ms1_id_end, 
                     FeatureParaPtr para_ptr) {
  if (ms1_ptr_vec.size() == 0) return;
  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  std::vector<double> ext_masses = para_ptr->getExtendMasses(prec_mass);
  for (int i = ms1_id_begin; i <= ms1_id_end; i++) {
    DeconvMsPtr ms1_ptr = ms1_ptr_vec[i];
    for (size_t j = 0; j < ms1_ptr->size(); j++) {
      DeconvPeakPtr peak = ms1_ptr->getPeakPtr(j);
      if (peak != nullptr) {
        for (size_t k = 0; k < ext_masses.size(); k++) {
          double mass_diff = std::abs(ext_masses[k] - peak->getPosition());
          if (mass_diff <= error_tole) {
            matched_peaks.push_back(peak);
            break;
          }
        }
      }
    }
  }
}

bool containPrecursor(DeconvMsPtr ms1_ptr, double prec_mass, FeatureParaPtr para_ptr) {
  if (ms1_ptr == nullptr) return false;
  // double prec_chrg = best_ptr->getPrecCharge();
  std::vector<double> ext_masses = para_ptr->getExtendMasses(prec_mass);

  double min_diff = std::numeric_limits<double>::max();
  for (size_t i = 0; i < ms1_ptr->size(); i++) {
    DeconvPeakPtr peak = ms1_ptr->getPeakPtr(i);
    // if (peak->getCharge() == prec_chrg) {
    // do not test charge
    if (peak != nullptr) {
      for (size_t j = 0; j < ext_masses.size(); j++) {
        double mass_diff = std::abs(ext_masses[j] - peak->getPosition());
        if (mass_diff < min_diff) {
          min_diff = mass_diff;
        }
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

int getMs1IdBegin(const DeconvMsPtrVec &ms1_ptr_vec, int sp_id, double prec_mass,
                  FeatureParaPtr para_ptr) {
  int cur_id = sp_id;
  int result_id = sp_id;
  int miss_num = 0;
  while (cur_id > 0) {
    if (containPrecursor(ms1_ptr_vec[cur_id], prec_mass, para_ptr)) {
      miss_num = 0;
      result_id = cur_id;
    } else {
      miss_num++;
    }
    cur_id--;
    if (miss_num >= 15) {
      break;
    }
  }
  return result_id;
}

int getMs1IdEnd(const DeconvMsPtrVec &ms1_ptr_vec, int sp_id, double prec_mass,
                FeatureParaPtr para_ptr) {
  int cur_id = sp_id;
  int result_id = sp_id;
  int miss_num = 0;
  while (cur_id < static_cast<int>(ms1_ptr_vec.size()) ) {
    if (containPrecursor(ms1_ptr_vec[cur_id], prec_mass, para_ptr)) {
      miss_num = 0;
      result_id = cur_id;
    } else {
      miss_num++;
    }
    cur_id++;
    if (miss_num >= 15) {
      break;
    }
  }
  return result_id;
}

double getFeatureInte(DeconvPeakPtrVec &matched_peaks) {
  double inte = 0;
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    inte += matched_peaks[i]->getIntensity();
  }
  return inte;
}

double getFeatureMass(double best_peak_mass, DeconvPeakPtrVec &matched_peaks,
                      FeatureParaPtr para_ptr) {
  std::vector<double> offsets = para_ptr->getExtendOffsets();
  size_t offset_num = offsets.size();
  DeconvPeakPtrVec2D offset_peaks;
  for (size_t i = 0; i < offsets.size(); i++) {
    DeconvPeakPtrVec peaks;
    offset_peaks.push_back(peaks);
  }
  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(best_peak_mass);
  // get peaks for each offset
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    DeconvPeakPtr peak = matched_peaks[i];
    for (size_t j = 0; j < offset_num; j++) {
      double mass_diff = std::abs(best_peak_mass + offsets[j] - peak->getPosition());
      if (mass_diff <= error_tole) {
        offset_peaks[j].push_back(peak);
        break;
      }
    }
  }
  // get feature intensities for each offset
  std::vector<double> intes(offset_num, 0);
  for (size_t i = 0; i < offset_num; i++) {
    for (size_t j = 0; j < offset_peaks[i].size(); j++) {
      intes[i] += offset_peaks[i][j]->getIntensity();
    }
  }
  // get the best offset
  int best_offset = -1;
  double best_inte = -1;
  for (size_t i = 0; i < offset_num; i++) {
    if (intes[i] > best_inte) {
      best_inte = intes[i];
      best_offset = i;
    }
  }
  //get the mass as weighted average
  double inte_sum = 0;
  double product_sum = 0;
  for (size_t i = 0; i < offset_peaks[best_offset].size(); i++) {
    DeconvPeakPtr peak = offset_peaks[best_offset][i];
    inte_sum += peak->getIntensity();
    product_sum = product_sum + peak->getPosition() * peak->getIntensity();
  }
  if (inte_sum == 0.0) {
    LOG_ERROR("Intensity sum is 0!");
  }
  return product_sum/inte_sum;
}

int getMinCharge (DeconvPeakPtrVec &matched_peaks) {
  int min_charge = std::numeric_limits<int>::max();
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    if (matched_peaks[i]->getCharge() < min_charge) {
      min_charge = matched_peaks[i]->getCharge();
    }
  }
  return min_charge;
}

int getMaxCharge (DeconvPeakPtrVec &matched_peaks) {
  int max_charge = -1;
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    if (matched_peaks[i]->getCharge() > max_charge) {
      max_charge = matched_peaks[i]->getCharge();
    }
  }
  return max_charge;
}

FracFeaturePtr getFeature(int sp_id, double prec_mass, int feat_id, 
                          DeconvMsPtrVec &ms1_ptr_vec,
                          DeconvPeakPtrVec &matched_peaks, FeatureParaPtr para_ptr) {
  int ms1_id_begin = getMs1IdBegin(ms1_ptr_vec, sp_id, prec_mass, para_ptr);
  int ms1_id_end = getMs1IdEnd(ms1_ptr_vec, sp_id, prec_mass, para_ptr);
  getMatchedPeaks(ms1_ptr_vec, prec_mass, matched_peaks,
                  ms1_id_begin, ms1_id_end, para_ptr);
  if (matched_peaks.size() == 0) {
    return nullptr;
  }
  double feat_inte = getFeatureInte(matched_peaks);
  double feat_mass = getFeatureMass(prec_mass, matched_peaks, para_ptr);
  int min_charge = getMinCharge(matched_peaks);
  int max_charge = getMaxCharge(matched_peaks);
  double retent_begin = ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr()->getRetentionTime();
  double retent_end = ms1_ptr_vec[ms1_id_end]->getMsHeaderPtr()->getRetentionTime();
  int ms1_scan_begin = ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr()->getFirstScanNum();
  int ms1_scan_end = ms1_ptr_vec[ms1_id_end]->getMsHeaderPtr()->getFirstScanNum();
  FracFeaturePtr feature_ptr = std::make_shared<FracFeature>(feat_id, 
                                                             para_ptr->frac_id_,
                                                             para_ptr->file_name_,
                                                             feat_mass,
                                                             feat_inte,
                                                             retent_begin,
                                                             retent_end,
                                                             ms1_scan_begin,
                                                             ms1_scan_end, min_charge,
                                                             max_charge);
  return feature_ptr;
}

void removePeaks (DeconvMsPtrVec &ms1_ptr_vec, DeconvPeakPtrVec &matched_peaks) {
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    DeconvPeakPtr peak = matched_peaks[i];
    int sp_id = peak->getSpId();
    int peak_id = peak->getId();
    ms1_ptr_vec[sp_id]->setPeakPtrNull(peak_id);
  }
}

void findMsOneFeatures(DeconvMsPtrVec &ms1_ptr_vec, FeatureParaPtr para_ptr,
                       FracFeaturePtrVec &features) {
  //get all peaks
  DeconvPeakPtrVec all_peaks;
  for (size_t i = 0; i < ms1_ptr_vec.size(); i++) {
    DeconvPeakPtrVec peaks = ms1_ptr_vec[i]->getPeakPtrVec();
    all_peaks.insert(std::end(all_peaks), std::begin(peaks), std::end(peaks));
  }
  //sort all peaks by intensities
  std::sort(all_peaks.begin(), all_peaks.end(), Peak::cmpInteDec);
  int feat_id = 0;
  size_t peak_idx = 0;
  while (static_cast<int>(features.size()) < para_ptr->feature_num_ 
         && peak_idx < all_peaks.size()) {
    DeconvPeakPtr best_peak = all_peaks[peak_idx];
    if (peakExists(ms1_ptr_vec, best_peak)) {
      //std::cout << "Find feature " << feat_id << " peak intensity " << best_peak->getIntensity() << std::endl; 
      int sp_id = best_peak->getSpId();
      double prec_mass = best_peak->getPosition();
      DeconvPeakPtrVec matched_peaks;
      FracFeaturePtr feature_ptr = getFeature(sp_id, prec_mass, feat_id, ms1_ptr_vec,
                                              matched_peaks, para_ptr);
      if (feature_ptr != nullptr) {
        features.push_back(feature_ptr);
        removePeaks(ms1_ptr_vec, matched_peaks);
        feat_id++;
      }
    }
    peak_idx++;
  }
}

void readHeaders(const std::string & file_name, MsHeaderPtrVec &header_ptr_vec) {
  int sp_num_in_group = 1;
  MsAlignReader sp_reader(file_name, sp_num_in_group, nullptr,
                          std::set<std::string>());
  DeconvMsPtr ms_ptr;
  LOG_DEBUG("Start search");
  while ((ms_ptr = sp_reader.getNextMs()) != nullptr) {
    header_ptr_vec.push_back(ms_ptr->getMsHeaderPtr());
    //std::cout << std::flush <<  "reading spectrum " << header_ptr_vec.size() << "\r";
  }
  sp_reader.close();
  //std::cout << std::endl;
}

inline bool isMatch(FracFeaturePtr feature_ptr, MsHeaderPtr header, FeatureParaPtr para_ptr) {
  int ms1_scan = header->getMsOneScan();
  if (ms1_scan < feature_ptr->getScanBegin()) {
    return false;
  }
  if (ms1_scan > feature_ptr->getScanEnd()) {
    return false;
  }
  double prec_mass = header->getPrecMonoMass();
  std::vector<double> search_masses = para_ptr->getSearchMasses(prec_mass);

  double feature_mass = feature_ptr->getMonoMass();

  double min_diff = std::numeric_limits<double>::max();
  for (size_t j = 0; j < search_masses.size(); j++) {
    double mass_diff = std::abs(search_masses[j] - feature_mass);
    if (mass_diff < min_diff) {
      min_diff = mass_diff;
    }
  }
  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  if (min_diff <= error_tole) {
    return true;
  }
  return false;
}

inline FracFeaturePtr getMatchedFeaturePtr(FracFeaturePtrVec &features, 
                                           MsHeaderPtr header,
                                           FeatureParaPtr para_ptr) {
  for (size_t i = 0; i < features.size(); i++) {
    if (isMatch(features[i], header, para_ptr)) {
      return features[i];
    }
  }
  return nullptr;
}

void getMs2Features(DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtrVec &header_ptr_vec, 
                    FracFeaturePtrVec &features, FeatureParaPtr para_ptr, 
                    SpecFeaturePtrVec &ms2_features) {
  MsHeaderPtrVec sorted_ptrs = header_ptr_vec;
  std::sort(sorted_ptrs.begin(), sorted_ptrs.end(), MsHeader::cmpPrecInteDec);
  for (size_t i = 0; i < sorted_ptrs.size(); i++) {
    MsHeaderPtr header = sorted_ptrs[i];
    FracFeaturePtr ft_ptr = getMatchedFeaturePtr(features, header, para_ptr);
    if (ft_ptr != nullptr) {
      SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(header, ft_ptr);
      ms2_features.push_back(ms2_feature);
    }
    else {
      int sp_id = header->getMsOneId();
      double prec_mass = header->getPrecMonoMass();
      if (prec_mass > 0) {
        int feat_id = static_cast<int>(features.size());
        DeconvPeakPtrVec matched_peaks;
        FracFeaturePtr feature_ptr = getFeature(sp_id, prec_mass, feat_id, ms1_ptr_vec,
                                                matched_peaks, para_ptr);
        // if we find a feature in ms1.msalign
        // it is possible that some ms headers do not have matched features. 
        if (feature_ptr != nullptr) {
          features.push_back(feature_ptr);
          removePeaks(ms1_ptr_vec, matched_peaks);

          SpecFeaturePtr ms2_feature = std::make_shared<SpecFeature>(header, feature_ptr);
          ms2_features.push_back(ms2_feature);
        }
        else {
          LOG_WARN("Cannot find features in LC/MS! Spectrum id: " << sp_id);
        }
      }
    }
  }
}

/*
void writeMs2Feature(const std::string & output_file_name,
                     const MsHeaderPtrVec &ms2_header_ptrs,
                     std::string argu_str) {
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
      << "FRACTION_FEATURE_ID" << "\t"
      << "FRACTION_FEATURE_INTENSITY" << std::endl;
  for (size_t i = 0; i < ms2_header_ptrs.size(); i++) {
    MsHeaderPtr header = ms2_header_ptrs[i];
    of << header->getId() << "\t"
        << header->getScansString() << "\t"
        << header->getMsOneId() << "\t"
        << header->getMsOneScan() << "\t"
        << header->getPrecMonoMass() << "\t"
        << header->getPrecInte() << "\t";
        */
    /*
    if (header->getFeatureId() >= 0) {
      of << header->getFeatureId() << "\t"
         << header->getFeatureInte() << std::endl;
    } 
    else {
      of << "-1" << "\t" << "-1" << std::endl;
    }
    */
/*
  }
  of.close();
}
*/

void process(int frac_id, std::string &sp_file_name, 
             bool missing_level_one, std::string &argu_str) {
  //logger::setLogLevel(2);
  FeatureParaPtr para_ptr = std::make_shared<FeaturePara>(frac_id, sp_file_name);
  std::string base_name = file_util::basename(sp_file_name);
  // read ms1 deconvoluted spectra
  std::string ms1_file_name = base_name + "_ms1.msalign";
  DeconvMsPtrVec ms1_ptr_vec;
  FracFeaturePtrVec features;
  if (!missing_level_one) {
    MsAlignReader::readMsOneSpectra(ms1_file_name, ms1_ptr_vec);
    findMsOneFeatures(ms1_ptr_vec, para_ptr, features);
  }

  LOG_DEBUG("start reading ms2");
  std::string ms2_file_name = base_name + "_ms2.msalign";
  MsHeaderPtrVec header_ptr_vec;
  readHeaders(ms2_file_name, header_ptr_vec);
  SpecFeaturePtrVec ms2_features;
  getMs2Features(ms1_ptr_vec, header_ptr_vec, features, para_ptr, ms2_features);

  std::sort(features.begin(), features.end(), FracFeature::cmpMassInc);
  std::string output_file_name = base_name + "_frac.feature";
  frac_feature_writer::writeFeatures(output_file_name, features);

  output_file_name = base_name + "_spec.feature";
  spec_feature_writer::writeFeatures(output_file_name, ms2_features); 

}

}  // namespace 

}  // namespace toppic 

