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
#include <algorithm>
#include "common/util/file_util.hpp"
#include "common/base/mod_util.hpp"
#include "seq/fasta_util.hpp"
#include "seq/fasta_index_reader.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_reader.hpp"
#include "deconv/feature/feature.hpp"
#include "deconv/feature/feature_para.hpp"
#include "quant/topquant_process.hpp"

namespace toppic {

namespace topquant_process {

void readSpectra(const std::string &file_name, DeconvMsPtrVec &ms_ptr_vec) {
  int sp_num_in_group = 1;
  MsAlignReader sp_reader(file_name, sp_num_in_group,
                          nullptr, std::set<std::string>());

  DeconvMsPtr ms_ptr;
  //LOG_DEBUG("Start search");
  while ((ms_ptr = sp_reader.getNextMs())!= nullptr) {
    ms_ptr->getMsHeaderPtr()->setMsLevel(1);
    ms_ptr_vec.push_back(ms_ptr);
    std::cout << std::flush <<  "reading spectrum " << ms_ptr_vec.size() << "\r";
  }
  sp_reader.close();
  std::cout << "reading spectrum finish" <<  std::endl;
}

bool containPrecursor(DeconvMsPtr ms1_ptr, double prec_mass, FeatureParaPtr para_ptr) {
  if (ms1_ptr == nullptr) return false;
  // double prec_chrg = best_ptr->getPrecCharge();
  std::vector<double> ext_masses = para_ptr->getExtMasses(prec_mass);

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

void getMatchedPeaks(DeconvMsPtrVec &ms1_ptr_vec, double prec_mass,
                     DeconvPeakPtrVec &matched_peaks, 
                     int ms1_id_begin, int ms1_id_end, 
                     FeatureParaPtr para_ptr) {
  if (ms1_ptr_vec.size() == 0) return;
  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
  std::vector<double> ext_masses = para_ptr->getExtMasses(prec_mass);
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

double getFeatureInte(DeconvPeakPtrVec &matched_peaks) {
  double inte = 0;
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    inte += matched_peaks[i]->getIntensity();
  }
  return inte;
}

double getFeatureMass(double best_peak_mass, DeconvPeakPtrVec &matched_peaks,
                      FeatureParaPtr para_ptr) {
  std::vector<double> offsets = para_ptr->getExtOffsets();
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

void removePeaks (DeconvMsPtrVec &ms1_ptr_vec, DeconvPeakPtrVec &matched_peaks) {
  for (size_t i = 0; i < matched_peaks.size(); i++) {
    DeconvPeakPtr peak = matched_peaks[i];
    int sp_id = peak->getSpId();
    int peak_id = peak->getId();
    ms1_ptr_vec[sp_id]->setPeakPtrNull(peak_id);
  }
}

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

FeaturePtr getFeature(int sp_id, double prec_mass, int feat_id, DeconvMsPtrVec &ms1_ptr_vec,
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
  FeaturePtr feature_ptr = std::make_shared<Feature>(feat_id, feat_mass,
                                                     feat_inte,
                                                     retent_begin,
                                                     retent_end,
                                                     ms1_scan_begin,
                                                     ms1_scan_end, min_charge,
                                                     max_charge);
  return feature_ptr;
}

int findBestFeature(FeaturePtrVec &features, PrsmStrPtr prsm, FeatureParaPtr para_ptr) {
  double mass = prsm->getOriPrecMass();
  double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(mass);
  int spec_scan = prsm->getSpectrumScan(); 
  int best_feature_id = -1;
  double best_error = std::numeric_limits<double>::max();
  std::vector<double> ext_offsets = para_ptr->getExtOffsets();
  for (size_t i = 0; i < features.size(); i++) {
    int feat_scan_begin = features[i]->getScanBegin();
    int feat_scan_end = features[i]->getScanEnd();
    if (spec_scan < feat_scan_begin || spec_scan > feat_scan_end + 10) {
      continue;
    }
    for (size_t k = 0; k < ext_offsets.size(); k++) {
      double mass_diff = std::abs(mass + ext_offsets[k] - features[i]->getMonoMass());
      if (mass_diff <= error_tole) {
        double cur_error = mass_diff + std::abs(ext_offsets[k]);
        if (cur_error < best_error) {
          best_feature_id = i;
          best_error = cur_error;
        }
      }
    }
  }
  return best_feature_id;
}

int findSpId(DeconvMsPtrVec &ms1_ptr_vec, int scan) {
  for (size_t i = 0; i < ms1_ptr_vec.size() -1; i++) {
    if (scan > ms1_ptr_vec[i]->getMsHeaderPtr()->getFirstScanNum()
        && scan < ms1_ptr_vec[i+1]->getMsHeaderPtr()->getFirstScanNum()) {
      return ms1_ptr_vec[i]->getMsHeaderPtr()->getId();
    }
  }
  return -1;
}

void findFeatures(DeconvMsPtrVec &ms1_ptr_vec, FeatureParaPtr para_ptr,
                  FeaturePtrVec &features, PrsmStrPtrVec &prsms) {

  for (size_t i = 0; i < prsms.size(); i++) {
    PrsmStrPtr prsm = prsms[i];
    int id = findBestFeature(features, prsm, para_ptr);
    // if not found 
    if (id <  0) {
      int sp_scan  = prsm->getSpectrumScan();
      int sp_id = findSpId(ms1_ptr_vec, sp_scan);
      if (sp_id < 0) {
        LOG_ERROR("Cannot find sp id!");
        continue;
      }
      double prec_mass = prsm->getOriPrecMass();
      DeconvPeakPtrVec matched_peaks;
      FeaturePtr feature_ptr = getFeature(sp_id, prec_mass, feat_id, ms1_ptr_vec,
                                          matched_peaks, para_ptr);
      if (feature_ptr != nullptr) {
        features.push_back(feature_ptr);
        removePeaks(ms1_ptr_vec, matched_peaks);
      }
      feat_id++;
    }
  }
  std::sort(features.begin(), features.end(), Feature::cmpMassInc);
}

void writeFeatures(const std::string &output_file_name,
                   const FeaturePtrVec &features, 
                   const PrsmStrPtrVec &prsms) {
  std::ofstream of(output_file_name, std::ofstream::out);
  of.precision(16);
  of << "ID" << "\t"
      << "Mass" << "\t"
      << "Intensity" << "\t"
      << "Time begin" << "\t"
      << "Time end" << "\t"
      << "First scan" << "\t"
      << "Last scan" << "\t"
      << "Minimum charge state" << "\t"
      << "Maximum charge state" << "\t"
      << "Protein" << "\t"
      << "Protein description " << "\t"
      << "First residue" << "\t" 
      << "Last residue" << "\t"
      << "Proteoform" << "\t"
      << "MS2 Scan" <<  "\t"
      << "Precursor mass" <<  std::endl;
  for (size_t i = 0; i < features.size(); i++) {
    FeaturePtr feature = features[i];
    of << feature->getId() << "\t"
        << feature->getMonoMass() << "\t"
        << feature->getIntensity() << "\t"
        << feature->getRetentBegin() << "\t"
        << feature->getRetentEnd() << "\t"
        << feature->getScanBegin() << "\t"
        << feature->getScanEnd() << "\t"
        << feature->getMinCharge() << "\t"
        << feature->getMaxCharge() << "\t";
    PrsmStrPtr prsm = prsms[i];
    if (prsm != nullptr) {
      of << prsm->getSeqName() << "\t"
          << prsm->getSeqDesc() << "\t"
          << (prsm->getProteoformStartPos() + 1) << "\t"
          << (prsm->getProteoformEndPos() + 1) << "\t"
          << prsm->getProteinMatchSeq() << "\t"
          << prsm->getSpectrumScan() << "\t"
          << prsm->getOriPrecMass() << std::endl;
    }
    else {
      of << "\t"
          << "\t"
          << "\t"
          << "\t"
          << "\t"
          << "\t"
          << std::endl;
    }
  }
  of.close();
}


PrsmStrPtrVec readPrsms(std::string &prsm_file_name) {
  std::string ori_db_file_name = "uniprot_zebrafish.fasta";
  std::string db_file_name = ori_db_file_name + "_target";
  fasta_util::dbSimplePreprocess(ori_db_file_name, db_file_name);
  std::string fixed_mod = "C57";

  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);
  ModPtrVec fix_mod_ptr_vec = mod_util::geneFixedModList(fixed_mod);
  PrsmStrPtrVec prsms = PrsmReader::readAllPrsmStrsMatchSeq(prsm_file_name,
                                                            seq_reader,
                                                            fix_mod_ptr_vec);
  return prsms;
}

void matchPrsms(FeaturePtrVec &features, PrsmStrPtrVec &matched_prsms, 
                PrsmStrPtrVec &prsms, FeatureParaPtr para_ptr) {
  for (size_t i = 0; i < prsms.size(); i++) {
    PrsmStrPtr prsm = prsms[i];
    int id = findBestFeature(features, prsm, para_ptr);
    if (id >= 0) {
      if (matched_prsms[id] == nullptr || matched_prsms[id]->getEValue() > prsm->getEValue()) {
        matched_prsms[id] = prsm;
      }
    }
    else {
      LOG_ERROR("Prsm not matched: mass " << prsm->getAdjustedPrecMass() << " scan " 
                << prsm->getSpectrumScan() << " inte " << prsm->getPrecFeatureInte());
    }
  }
}

void process(std::string &sp_file_name, bool missing_level_one, 
             std::string &argu_str) {
  //logger::setLogLevel(2);
  /*
  FeatureParaPtr para_ptr = std::make_shared<FeaturePara>();
  // read ms1 deconvoluted spectra
  std::string ms1_file_name = base_name + "_ms1.msalign";
  DeconvMsPtrVec ms1_ptr_vec;
  if (!missing_level_one) readSpectra(ms1_file_name, ms1_ptr_vec);
  */
  FeaturePtrVec features;
  std::string base_name = file_util::basename(sp_file_name);
  std::string prsm_file_name = base_name + "_ms2_toppic_proteoform.xml";
  PrsmStrPtrVec prsms = readPrsms(prsm_file_name);
  std::cout << "start finding feature" << std::endl;
  findFeatures(ms1_ptr_vec, para_ptr, features, prsms);
  std::cout << "end finding feature" << std::endl;

  PrsmStrPtrVec matched_prsms(features.size());
  matchPrsms(features, matched_prsms, prsms, para_ptr);

  std::string output_file_name = base_name + ".feature";
  writeFeatures(output_file_name, features, matched_prsms);
}

}  // namespace feature_detect_2
*/

}  // namespace toppic 
