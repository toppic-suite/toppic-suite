//
// Created by abbash on 8/23/22.
//

#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "peak_matrix.hpp"
#include "ms/spec/peak.hpp"

toppic::PeakMatrix::PeakMatrix(PeakPtrVec2D& raw_peaks, DeconvMsPtrVec ms1_ptr_vec, double bin_size, double snr){
  /// get min mz value
  std::vector<double> intes;
  std::vector<double> mz;
  for (auto & raw_peak : raw_peaks){
    for (auto & p : raw_peak) {
      mz.push_back(p->getPosition());
      intes.push_back(p->getIntensity());
    }
  }
  /// set values!
  min_mz_ = *std::min_element(mz.begin(), mz.end());
  max_mz_ = *std::max_element(mz.begin(), mz.end());
//  min_inte_ = getDataLevelNoiseIntensities(intes);
  min_inte_ = 188.80874;
  specs_ = get_spec_list(ms1_ptr_vec);
  bin_size_ = bin_size;
  spec_num_ = raw_peaks.size();
  bin_num_ = int((max_mz_ - min_mz_)/bin_size_) + 1;
  std::cout << "Data Level Noise Intensity Level: " << min_inte_ << std::endl;
  std::cout << "Min and Max m/z values: " << min_mz_ << " , " << max_mz_ << std::endl;
  init_matrix(raw_peaks, snr);
}

toppic::spec_list toppic::PeakMatrix::get_spec_list(DeconvMsPtrVec ms1_ptr_vec){
  std::vector<Spectrum> spec_list;
  for (auto & i : ms1_ptr_vec)
    spec_list.emplace_back(i->getMsHeaderPtr()->getId(), i->getMsHeaderPtr()->getFirstScanNum(), i->getMsHeaderPtr()->getRetentionTime());
  return spec_list;
}

void toppic::PeakMatrix::init_matrix(PeakPtrVec2D& raw_peaks, double snr){
//  std::map<int, PeakRow> peak_matrix;
  for (int spec_id = 0; spec_id < spec_num_; spec_id++){
    if (spec_id%100 == 0)
      std::cout << "Init matrix processing peaks in spectrum " << spec_id << " out of " << spec_num_ << std::endl;
    // Fill the peak matrix

    matrix_[spec_id] = PeakRow(specs_[spec_id], bin_num_);

    std::vector<PeakPtr> spec_peak_list = raw_peaks[spec_id];
    std::vector<ExpPeak> exp_peak_list;
    for (auto& p : spec_peak_list){
      int peak_id = peaks_.size();
      ExpPeak new_peak = ExpPeak(peak_id, spec_id, p);
      exp_peak_list.push_back(new_peak);
      peaks_.push_back(new_peak);
    }
    for (auto& cur_peak : exp_peak_list) {
      // Only keep peaks above data level noise intensity * SNR
      if (cur_peak.getInte() > snr * min_inte_) {
        int bin_idx = get_index(cur_peak.getPos());
        matrix_[spec_id].addPeak(bin_idx, cur_peak);
      }
    }
//    matrix_ = peak_matrix;
    //std::cout << "peak matrix size " << matrix_.size() << " peak num " << peaks_.size() << std::endl;
  }
}

int toppic::PeakMatrix::get_index(double mz) {
  double mz_diff = mz - min_mz_;
  int bin_idx = int(mz_diff /bin_size_);
  return bin_idx;
}

void toppic::PeakMatrix::find_pair_neighbors(int spec_id, int search_bin_num, double mass_tol){
  std::vector<std::vector<ExpPeak>> first_bin_list = matrix_[spec_id].getRow();
  std::vector<std::vector<ExpPeak>> second_bin_list = matrix_[spec_id+1].getRow();
  for (int bin_idx = 0; bin_idx < bin_num_; bin_idx++) {
    int start = std::max(0, bin_idx - search_bin_num);
    int end = std::min(bin_idx + search_bin_num, bin_num_ - 1);
  //  std::cout << "bin_info: " << bin_idx << ", " << first_bin_list[bin_idx].size() << ", " << start << ", " << end << std::endl;
    for (auto & first_peak : first_bin_list[bin_idx]) {
      for (int second_idx = start; second_idx < end + 1; second_idx++) {
        for (auto & second_peak : second_bin_list[second_idx]) {
          double mass_diff = std::abs(first_peak.getPos() - second_peak.getPos());
          if (mass_diff <= mass_tol) {
  //          std::cout << "inside loop 3 " << first_peak.getPos() << " " << second_peak.getPos() << " " << mass_diff << " " << (mass_diff <= mass_tol) << std::endl;
            first_peak.setNeighbor(true);
            second_peak.setNeighbor(true);
//            break;
          }
        }
      }
    }
  }
  matrix_[spec_id].setRow(first_bin_list);
  matrix_[spec_id + 1].setRow(second_bin_list);
}

void toppic::PeakMatrix::find_remove_non_neighbors(double mass_tol) {
  int search_bin_num = int(mass_tol/bin_size_) + 1;
  for (auto & peak : peaks_)
    peak.setNeighbor(false);
  for (int spec_id = 0; spec_id < spec_num_ - 1; spec_id++)
    find_pair_neighbors(spec_id, search_bin_num, mass_tol);
  // remove peaks without neighbors
  for (int spec_id = 0; spec_id < spec_num_ - 1; spec_id++){
    std::vector<std::vector<ExpPeak>> bin_list = matrix_[spec_id].getRow();
    for (int bin_idx = 0; bin_idx < bin_num_; bin_idx++)
      for (auto &peak: bin_list[bin_idx])
        if (!peak.getNeighbor())
          peak = ExpPeak();
    matrix_[spec_id].setRow(bin_list);
  }
}

void toppic::PeakMatrix::count_peaks(){
  int counter = 0;
  for (int spec_id = 0; spec_id < spec_num_ - 1; spec_id++){
    std::vector<std::vector<ExpPeak>> bin_list = matrix_[spec_id].getRow();
    for (const auto& bin: bin_list)
      for (auto peak: bin)
        if (!peak.isEmpty())
          counter++;
  }
  std::cout << "Number of Peaks: " << counter << std::endl;
}

void toppic::PeakMatrix::remove_peak(ExpPeak& peak){
  int spec_id = peak.getSpecId();
  int bin_idx = get_index(peak.getPos());

  std::vector<std::vector<ExpPeak>> peaks_row = matrix_[spec_id].getRow();
  std::vector<ExpPeak> bin_peaks = peaks_row[bin_idx];

  std::vector<ExpPeak> new_bin;
  for (auto bin_peak : bin_peaks) {
    if (bin_peak.getPeakId() == peak.getPeakId()) continue;
    new_bin.push_back(bin_peak);
  }
  peaks_row[bin_idx] = new_bin;
  matrix_[spec_id].setRow(peaks_row);
}

//void toppic::PeakMatrix::remove_peak(ExpPeak& peak){
//  int spec_id = peak.getSpecId();
//  int bin_idx = get_index(peak.getPos());
//  std::vector<std::vector<ExpPeak>> peaks_row = matrix_[spec_id].getRow();
//  std::vector<ExpPeak> bin_peaks = peaks_row[bin_idx];
//  for (auto bin_peak = bin_peaks.begin(); bin_peak < bin_peaks.end(); bin_peak++) {
//    if (bin_peak->getPeakId() == peak.getPeakId()) {
//      bin_peaks.erase(bin_peak);
//      break;
//    }
//  }
//  matrix_[spec_id].setRow(peaks_row);
//}

void toppic::PeakMatrix::remove_peak_in_range(int spec_id, double min_pos, double max_pos) {
  int start_bin_idx = get_index(min_pos);
  int end_bin_idx = get_index(max_pos);
  std::vector<std::vector<ExpPeak>> peaks_row = matrix_[spec_id].getRow();
  for (int bin_idx = start_bin_idx; bin_idx < end_bin_idx + 1; bin_idx++) {
    std::vector<ExpPeak> bin_peaks = peaks_row[bin_idx];
    for (auto bin_peak = bin_peaks.begin(); bin_peak < bin_peaks.end(); bin_peak++) {
      if (bin_peak->getPos() >= min_pos and bin_peak->getPos() <= max_pos)
        bin_peaks.erase(bin_peak);
    }
  }
  matrix_[spec_id].setRow(peaks_row);
}
