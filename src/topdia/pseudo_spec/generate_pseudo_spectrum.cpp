// Copyright (c) 2014 - 2023, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#include <cstddef>

#include "topdia/pseudo_spec/pseudo_spectrum.hpp"
#include "topdia/pseudo_spec/generate_pseudo_spectrum.hpp"

namespace toppic {

GeneratePseudoSpectrum::GeneratePseudoSpectrum(TopfdParaPtr topfd_para_ptr, 
                                               TopdiaParaPtr topdia_para_ptr) {
  std::string output_base_name = topfd_para_ptr->getOutputBaseName();

  DeconvMsPtrVec deconv_ms1_ptr_vec;
  std::string ms1_file_name = output_base_name + "_ms1.msalign";
  //std::cout << "Reading ms1 file " << ms1_file_name << std::endl;
  msalign_reader_util::readAllSpectra(ms1_file_name, deconv_ms1_ptr_vec);
  for (const auto &ms1_data : deconv_ms1_ptr_vec) {
    rt_ms1_.push_back(ms1_data->getMsHeaderPtr()->getRetentionTime() / 60);
    ms1_scan_.push_back(ms1_data->getMsHeaderPtr()->getFirstScanNum());
  }

  DeconvMsPtrVec deconv_ms2_ptr_vec;
  std::string ms2_file_name = output_base_name + "_raw_ms2.msalign";
  //std::cout << "Reading ms2 file " << ms2_file_name << std::endl;
  msalign_reader_util::readAllSpectra(ms2_file_name, deconv_ms2_ptr_vec);

  // get isolation window base mz
  std::set<std::pair<double,double>> win_set;
  for (auto &ms2_data: deconv_ms2_ptr_vec) {
    if (ms2_data->getMsHeaderPtr()->getMsLevel() == 1)
      continue;
    std::pair<double, double> cur_win(ms2_data->getMsHeaderPtr()->getPrecWinBegin(), 
                                      ms2_data->getMsHeaderPtr()->getPrecWinEnd()); 
    win_set.insert(cur_win);
  }
  std::vector<std::pair<double,double>> win_list(win_set.begin(), win_set.end()); 
  win_list_ = win_list;

  // get retention times
  for (auto cur_win : win_list_) {
    std::vector<double> rt_ms2_window;
    ActivationPtrVec activation_ptr_vec;
    for (auto &ms2_data : deconv_ms2_ptr_vec) {
      double ms_win_begin = ms2_data->getMsHeaderPtr()->getPrecWinBegin();
      if (ms_win_begin != cur_win.first) continue;
      rt_ms2_window.push_back(ms2_data->getMsHeaderPtr()->getRetentionTime() / 60);
      activation_ptr_vec.push_back(ms2_data->getMsHeaderPtr()->getActivationPtr());
    }
    while (rt_ms2_window.size() < rt_ms1_.size()) {
      rt_ms2_window.push_back(rt_ms2_window.back());
      activation_ptr_vec.push_back(activation_ptr_vec.back());
    }
    rt_ms2_.push_back(rt_ms2_window);
    activation_ms2_.push_back(activation_ptr_vec);
  }

  // read feature files
  std::string filename = output_base_name + "_frac_ms1.mzrt.csv";
  //std::cout << "Reading ms1 feature file " << filename << std::endl;
  ms1_features_ = MzrtFeature::read_record(filename);
  for (auto cur_win : win_list_) {
    filename = output_base_name + "_" + std::to_string(cur_win.first) 
               + "_frac_ms2.mzrt.csv";
    //std::cout << "Reading ms2 feature file " << filename << std::endl;
    MzrtFeaturePtrVec feature_list = MzrtFeature::read_record(filename);
    ms2_features_.push_back(feature_list);
  }

  // get interpolated xic
  double max_rt = get_max_rt();
  double interval = 0.01;
  std::vector<double> rt_target;
  double i = 0;
  while (i <= max_rt) {
    rt_target.push_back(i);
    i += interval;
  }

  
  for (auto &ms1_feature : ms1_features_)
    ms1_feature->setInterpolatedXic(
        interp(rt_target, rt_ms1_, ms1_feature->getXic()));

  for (std::size_t iso_win_idx = 0; iso_win_idx < win_list_.size(); iso_win_idx++) {
    for (auto &ms2_feature : ms2_features_[iso_win_idx]) {
      ms2_feature->setInterpolatedXic(
          interp(rt_target, rt_ms2_[iso_win_idx], ms2_feature->getXic()));
    }
  }
}

void GeneratePseudoSpectrum::process(TopfdParaPtr topfd_para_ptr, 
                                     TopdiaParaPtr topdia_para_ptr) {
  int feature_id = 0;
  EnvParaPtr env_para_ptr =
      std::make_shared<EnvPara>(topfd_para_ptr->getMzError());
  std::string output_base_name = topfd_para_ptr->getOutputBaseName();
  std::string ms2_msalign_name = output_base_name + "_ms2.msalign";
  std::ofstream output;
  output.open(ms2_msalign_name); 

  for (std::size_t iso_win_idx = 0; iso_win_idx < win_list_.size(); iso_win_idx++) {
    std::pair<double, double> win = win_list_[iso_win_idx];
    MzrtFeaturePtrVec selected_ms1_features = get_iso_win_ms1_features(iso_win_idx);
    std::sort(selected_ms1_features.begin(), selected_ms1_features.end(),
              compareFeaturesInte);
    std::cout << "Processing isolation window [" << win.first << "," << win.second << "] with "
              << selected_ms1_features.size() << " features." << std::endl;
    for (auto &ms1_feature : selected_ms1_features) {
      int apex_cycle_distance_tole =
          std::min(3, ms1_feature->getCycleSpan() / 2);
      //TO CONTINUE
      ms1_feature->setWin(win);
      std::vector<PseudoPeak> pseudo_peak_list;
      for (std::size_t ms2_feature_idx = 0; ms2_feature_idx < ms2_features_[iso_win_idx].size();
           ms2_feature_idx++) {
        if (ms2_features_[iso_win_idx][ms2_feature_idx]->getUsedStatus())
          continue;
        MzrtFeaturePtr ms2_feature = ms2_features_[iso_win_idx][ms2_feature_idx];
        if (ms2_feature->getMass() < ms1_feature->getMass()) {
          int apex_cycle_distance =
              get_apex_cycle_distance(ms1_feature, ms2_feature);
          if (apex_cycle_distance > apex_cycle_distance_tole) continue;

          /// Get shared intensity
          std::vector<double> feature_xic_ms1 = ms1_feature->getInterpolatedXic();
          std::vector<double> feature_xic_ms2 = ms2_feature->getInterpolatedXic();
          double shared_area =
              computeSharedArea(MzrtFeature::normalizeXIC(feature_xic_ms1),
                                MzrtFeature::normalizeXIC(feature_xic_ms2));

          PseudoPeak peak = PseudoPeak(
              ms2_feature->getMass(), ms2_feature->getMonoMz(),
              ms2_feature->getCharge(), ms2_feature->getIntensity(), 0, 0,
              shared_area, ms2_feature->getCycleSpan(), apex_cycle_distance,
              ms2_feature->getTimeBegin(), ms2_feature->getTimeEnd(),
              ms2_feature->getApexCycle());
          peak.setMs2FeatureIdx(ms2_feature_idx);
          pseudo_peak_list.push_back(peak);
        }
      }
      /// filter peaks based on low-high mass dividers
      score_pseudo_peaks(pseudo_peak_list, ms1_feature);
      std::vector<PseudoPeak> filtered_pseudo_peak_list = filterPseudoPeaks(
          env_para_ptr, ms1_feature, ms2_features_[iso_win_idx],
          pseudo_peak_list, topdia_para_ptr->getPseudoScoreCutoff(),
          topdia_para_ptr->getPseudoMinPeaks());
      ms1_feature->setPseudoPeakNum(
          static_cast<int>(filtered_pseudo_peak_list.size()));
      PseudoSpectrumPtr pseudo_spec_ptr = std::make_shared<PseudoSpectrum>(
          ms1_feature, filtered_pseudo_peak_list);
      writePseudoSpectrum(output, topfd_para_ptr, topdia_para_ptr, 
                          feature_id, pseudo_spec_ptr->getMs1Feature(),
                          pseudo_spec_ptr->getFragmentFeatures(), iso_win_idx);
      feature_id++;
    }
  }
  output.close();
}

double GeneratePseudoSpectrum::get_max_rt() {
  double max_rt = rt_ms1_.back();
  for (int iso_win_idx = 0; iso_win_idx < win_list_.size(); iso_win_idx++) {
    if (max_rt > rt_ms2_[iso_win_idx].back())
      max_rt = rt_ms2_[iso_win_idx].back();
  }
  return max_rt;
}

std::vector<double> GeneratePseudoSpectrum::interp(
    const std::vector<double> &x, const std::vector<double> &xp,
    const std::vector<double> &fp) {
  std::vector<double> interpolatedValues;
  interpolatedValues.reserve(x.size());

  for (double xi : x) {
    auto it = std::lower_bound(xp.begin(), xp.end(), xi);
    if (it == xp.begin() || it == xp.end()) {
      interpolatedValues.push_back(0.0);
    } else {
      std::size_t i = it - xp.begin();
      double x1 = xp[i - 1];
      double x2 = xp[i];
      double y1 = fp[i - 1];
      double y2 = fp[i];
      double y = y1 + ((y2 - y1) / (x2 - x1)) * (xi - x1);
      interpolatedValues.push_back(y);
    }
  }
  return interpolatedValues;
}

MzrtFeaturePtrVec GeneratePseudoSpectrum::get_iso_win_ms1_features(int isolation_window_base_index) {
  MzrtFeaturePtrVec selected_features;
  std::pair<double, double> win = win_list_[isolation_window_base_index];
  double window_start_mz = win.first; 
  double window_end_mz = win.second;
  for (const auto &ms1_feature : ms1_features_) {
    std::vector<double> envelope_intensity = ms1_feature->getEnvelopeInte();
    std::vector<double> envelope_mz = ms1_feature->getEnvelopeMz();
    std::vector<double> xic = ms1_feature->getXic();
    double envelope_intensity_sum = std::accumulate(
        envelope_intensity.begin(), envelope_intensity.end(), 0.0);
    if (envelope_intensity_sum == 0.0) {
      continue;
    }
    double xic_sum = std::accumulate(xic.begin(), xic.end(), 0.0); 
    if (xic_sum == 0.0) {
      continue;
    }
  
    // get envelope peaks intensity present within isolation window
    double env_inte_sum_iso_win = 0;
    for (std::size_t idx = 0; idx < envelope_intensity.size(); idx++) {
      if (window_start_mz <= envelope_mz[idx] &&
          envelope_mz[idx] <= window_end_mz)
        env_inte_sum_iso_win += envelope_intensity[idx];
    }
    double coverage = env_inte_sum_iso_win / envelope_intensity_sum;
    if (coverage > 0.5) {
      selected_features.push_back(ms1_feature);
    }
  }
  std::sort(selected_features.begin(), selected_features.end(),
            compareFeaturesInte);
  return selected_features;
}

int GeneratePseudoSpectrum::get_apex_cycle_distance(
    const MzrtFeaturePtr &ms1_feature, const MzrtFeaturePtr &ms2_feature) {
  int apex_cycle_distance_tole = std::min(3, ms1_feature->getCycleSpan() / 2);
  int apex_cycle_distance =
      std::abs(ms1_feature->getApexCycle() - ms2_feature->getApexCycle());
  if (apex_cycle_distance > apex_cycle_distance_tole) {
    std::vector<double> ms1_xic = ms1_feature->getXic();
    std::vector<double> ms2_xic = ms2_feature->getXic();
    if (ms2_feature->getCycleSpan() > 15)
      ms2_xic = moving_avg(ms2_feature->getXic(), 3);

    int ms1_apex_scan = static_cast<int>(std::distance(
        ms1_xic.begin(), std::max_element(ms1_xic.begin(), ms1_xic.end())));
    int ms2_apex_scan = static_cast<int>(std::distance(
        ms2_xic.begin(), std::max_element(ms2_xic.begin(), ms2_xic.end())));
    apex_cycle_distance = std::abs(ms1_apex_scan - ms2_apex_scan);
    return apex_cycle_distance;
  }
  return apex_cycle_distance;
}

std::vector<double> GeneratePseudoSpectrum::moving_avg(std::vector<double> xic,
                                                       int size) {
  std::vector<double> smoothed_inte_list;
  std::vector<double> left_padding(1, 0);
  xic.insert(xic.begin(), left_padding.begin(), left_padding.end());
  int num_spec = static_cast<int>(xic.size());
  double sum = 0.0;
  int cnt = 0;
  for (int i = 0; i < num_spec; i++) {
    sum += xic[i];
    cnt++;
    if (cnt >= size) {
      smoothed_inte_list.push_back((sum / (double)size));
      sum -= xic[cnt - size];
    }
  }
  return smoothed_inte_list;
}

double GeneratePseudoSpectrum::computeSharedArea(
    const std::vector<double> &xic1, const std::vector<double> &xic2) {
  double sharedArea = 0.0;
  for (std::size_t i = 0; i < xic1.size(); ++i) {
    sharedArea += std::min(xic1[i], xic2[i]);
  }
  return sharedArea;
}

void GeneratePseudoSpectrum::score_pseudo_peaks(
    std::vector<PseudoPeak> &pseudo_peak_list,
    const MzrtFeaturePtr &ms1_feature) {
  std::sort(pseudo_peak_list.begin(), pseudo_peak_list.end(),
            comparePseudoPeaksInte);
  int total_peaks = static_cast<int>(pseudo_peak_list.size());
  double rank = total_peaks;
  for (int peak_idx = 0; peak_idx < total_peaks; peak_idx++) {
    pseudo_peak_list[peak_idx].setRank(rank / total_peaks);
    double score =
        get_pred(rank / total_peaks, pseudo_peak_list[peak_idx].getSharedInte(),
                 (pseudo_peak_list[peak_idx].getMS2CycleSpan() * 1.0) /
                     (ms1_feature->getCycleSpan() * 1.0));
    pseudo_peak_list[peak_idx].setScore(score);
    rank--;
  }
}

double GeneratePseudoSpectrum::get_pred(double intensity_ratio,
                                        double shared_area,
                                        double length_ratio) {
  double B0 = -3.349924626238689;
  double B1 = 1.8679961011204878;
  double B2 = 0.27006086659334383;
  double B3 = 3.983766800414337;
  double y = B0 + B1 * intensity_ratio + B2 * (length_ratio) + B3 * shared_area;

  double py = 1 / (1 + std::exp(-y));
  return py;
}

std::vector<PseudoPeak> GeneratePseudoSpectrum::filterPseudoPeaks(
    const EnvParaPtr &env_para_ptr, const MzrtFeaturePtr &ms1_feature,
    MzrtFeaturePtrVec &ms2_features_window,
    std::vector<PseudoPeak> &pseudo_peak_list, double cutoff,
    int min_peak_num) {
  std::sort(pseudo_peak_list.begin(), pseudo_peak_list.end(),
            comparePseudoPeaksScore);
  std::vector<PseudoPeak> low_mass_features;
  std::vector<PseudoPeak> high_mass_features;

  int low_mass_num = env_para_ptr->compLowMassNum();
  int high_mass_num = env_para_ptr->compHighMassNum(ms1_feature->getMass());

  // add shorter features with correlation
  int counter = 0;
  for (const auto &i : pseudo_peak_list) {
    if (i.getScore() < cutoff and counter >= min_peak_num) continue;
    if (i.getMass() <= env_para_ptr->low_high_dividor_) {
      if ((int)low_mass_features.size() < low_mass_num) {
        low_mass_features.push_back(i);
        ms2_features_window[i.getMs2FeatureIdx()]->setUsedStatus(true);
        counter++;
      }
    } else if ((int)high_mass_features.size() < high_mass_num) {
      high_mass_features.push_back(i);
      ms2_features_window[i.getMs2FeatureIdx()]->setUsedStatus(true);
      counter++;
    }
  }

  std::vector<PseudoPeak> result;
  result.insert(std::end(result), std::begin(low_mass_features),
                std::end(low_mass_features));
  result.insert(std::end(result), std::begin(high_mass_features),
                std::end(high_mass_features));
  std::sort(result.begin(), result.end(), comparePseudoPeaksScore);
  return result;
}

void GeneratePseudoSpectrum::writePseudoSpectrum(
    std::ofstream &output, TopfdParaPtr topfd_para_ptr, 
    TopdiaParaPtr topdia_para_ptr, 
    int ms1_feature_idx, MzrtFeaturePtr ms1_feature,
    std::vector<PseudoPeak> &assigned_ms2_features, 
    int iso_win_idx) {

  int ms1_apex_cycle = ms1_feature->getApexCycle();
  output << std::fixed;

  output << "BEGIN IONS" << std::endl;
  output << "FILE_NAME=" << topfd_para_ptr->getMzmlFileName() << std::endl;
  output << "SPECTRUM_ID=" << ms1_feature_idx << std::endl;  ///
  output << "TITLE=" << "Pseudo_Scan_" + std::to_string(ms1_feature_idx)
         << std::endl;
  output << "SCANS=" << (ms1_scan_[ms1_apex_cycle]+ 1 + iso_win_idx) << std::endl;  ///
  output << "RETENTION_TIME=" << std::fixed << std::setprecision(2)
         << ms1_feature->getTimeApex() << std::endl;
  output << "LEVEL=" << 2 << std::endl;
  output << "MS_ONE_ID=" << ms1_apex_cycle << std::endl;        ///
  // need to add MS_ONE_SCAN information                                                                          
  output << "MS_ONE_SCAN=" << ms1_scan_[ms1_apex_cycle] << std::endl;  ///
  output << "PRECURSOR_WINDOW_BEGIN=" << ms1_feature->getWin().first << std::endl;
  output << "PRECURSOR_WINDOW_END=" << ms1_feature->getWin().second << std::endl;
  // need to add ACTIVATION information                                                                          
  output << "ACTIVATION=" << activation_ms2_[iso_win_idx][ms1_apex_cycle]->getName() << std::endl;
  output << "PRECURSOR_MZ=" << std::setprecision(5) << ms1_feature->getMonoMz() << std::endl;
  output << "PRECURSOR_CHARGE=" << ms1_feature->getCharge() << std::endl;
  output << "PRECURSOR_MASS=" << ms1_feature->getMass() << std::endl;
  output << "PRECURSOR_INTENSITY=" << std::setprecision(2) << ms1_feature->getIntensity() << std::endl;
  // need to add MS1 feature ID                                                                           
  output << "PRECURSOR_FEATURE_ID=" << ms1_feature->getId() << std::endl;  ///
  output << "PRECURSOR_LENGTH=" << ms1_feature->getCycleSpan() << std::endl;

  for (const auto &peak : assigned_ms2_features) {
    output << std::fixed << std::setprecision(5) << peak.getMass();
    output << "\t" << std::fixed << std::setprecision(2) << peak.getIntensity();
    output << "\t" << peak.getCharge();
    output << "\t" << std::fixed << std::setprecision(2) << peak.getScore();
    output << "\t" << std::fixed << std::setprecision(2)
           << peak.getApexDiffScan();
    output << "\t" << std::fixed << std::setprecision(2) << peak.getRank();
    output << "\t" << std::fixed << std::setprecision(2)
           << peak.getMS2CycleSpan();
    output << "\t" << std::fixed << std::setprecision(2)
           << peak.getSharedInte();
    output << "\t" << peak.getMs2ApexCycle();
    output << std::endl;
  }
  output << "END IONS" << std::endl;
  output << std::endl;
}
}  // namespace toppic
