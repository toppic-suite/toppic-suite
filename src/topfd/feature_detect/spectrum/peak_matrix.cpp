//
// Created by abbash on 8/23/22.
//

#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "peak_matrix.hpp"
#include "ms/spec/peak.hpp"

toppic::PeakMatrix::PeakMatrix(PeakPtrVec2D raw_peaks, DeconvMsPtrVec ms1_ptr_vec){
  /// Input params
  double bin_size = 0.1;
  double snr = 3.0;


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
  min_inte_ = getDataLevelNoiseIntensities(intes);
  std::cout << "Data Level Noise Intensity Level: " << min_inte_ << std::endl;
  std::cout << "Min and Max m/z values: " << min_mz_ << " , " << max_mz_ << std::endl;

  specs_ = get_spec_list(ms1_ptr_vec);
  bin_size_ = 0.1;
  spec_num_ = raw_peaks.size();
  bin_num_ = int((max_mz_ - min_mz_)/bin_size_) + 1;
  init_matrix(raw_peaks, snr);
}

toppic::spec_list toppic::PeakMatrix::get_spec_list(DeconvMsPtrVec ms1_ptr_vec){
  std::vector<Spectrum> spec_list;
  for (int i = 0; i < ms1_ptr_vec.size(); i++){
    Spectrum spec = Spectrum(ms1_ptr_vec[i]->getMsHeaderPtr()->getMsOneId(), ms1_ptr_vec[i]->getMsHeaderPtr()->getMsOneScan(), ms1_ptr_vec[i]->getMsHeaderPtr()->getRetentionTime());
    spec_list.push_back(spec);
  }
  return spec_list;
}

void toppic::PeakMatrix::init_matrix(PeakPtrVec2D raw_peaks, double snr){
  std::map<int, PeakRow> peak_matrix;
  for (int spec_id = 0; spec_id < spec_num_; spec_id++){
    if (spec_id%100 == 0)
      std::cout << "Init matrix processing peaks in spectrum " << spec_id << " out of " << spec_num_ << std::endl;
    // Fill the peak matrix
    std::vector<PeakPtr> spec_peak_list = raw_peaks[spec_id];
    peak_matrix[spec_id] = PeakRow(specs_[spec_id], bin_num_);
    std::vector<ExpPeak> exp_peak_list;
    for (int p_idx = 0; p_idx < spec_peak_list.size(); p_idx++){
      PeakPtr p = spec_peak_list[p_idx];
      int peak_id = peaks_.size();
      ExpPeak new_peak = ExpPeak(peak_id, spec_id, p);
      exp_peak_list.push_back(new_peak);
      peaks_.push_back(new_peak);
    }
    for (int p_idx = 0; p_idx < exp_peak_list.size(); p_idx++) {
      ExpPeak cur_peak = exp_peak_list[p_idx];
      if (cur_peak.getInte() > snr * min_inte_) {
        // Only keep peaks above data level noise intensity * SNR
        int bin_idx = get_index(cur_peak.getPos());
        peak_matrix[spec_id].addPeak(bin_idx, cur_peak);
      }
    }
    matrix_ = peak_matrix;
    //std::cout << "peak matrix size " << matrix_.size() << " peak num " << peaks_.size() << std::endl;
  }
}

int toppic::PeakMatrix::get_index(double mz){
  double mz_diff = mz - min_mz_;
  int bin_idx = int(mz_diff /bin_size_);
  return bin_idx;
}

void toppic::PeakMatrix::find_pair_neighbors(PeakRow first_row, PeakRow second_row, int search_bin_num, double mass_tol){
  std::vector<std::vector<ExpPeak>> first_bin_list = first_row.getRow();
  std::vector<std::vector<ExpPeak>> second_bin_list = second_row.getRow();
  for (int bin_idx = 0; bin_idx < bin_num_; bin_idx++){
    int start = std::max(0, bin_idx - search_bin_num);
    int end = std::min(bin_idx + search_bin_num, bin_num_ - 1);
    for (int first_peak_idx = 0; first_peak_idx < first_bin_list[bin_idx].size(); first_peak_idx++) {
      ExpPeak first_peak = first_bin_list[bin_idx][first_peak_idx];
      for (int second_idx = start; second_idx < end + 1; second_idx++) {
        for (int second_peak_idx = 0; second_peak_idx < second_bin_list[bin_idx].size(); second_peak_idx++) {
          ExpPeak second_peak = second_bin_list[bin_idx][second_peak_idx];
          double mass_diff = abs(first_peak.getPos() - second_peak.getPos());
          if (mass_diff <= mass_tol) {
            first_peak.setNeighbor(true);
            second_peak.setNeighbor(true);
          }
        }
      }
    }
  }
}

void toppic::PeakMatrix::find_all_neighbors(double mass_tol){
  int search_bin_num = int(mass_tol/bin_size_) + 1;
  for (int i = 0; i < peaks_.size(); i++){
    ExpPeak peak = peaks_[i];
    peak.setNeighbor(false);
  }
  for (int spec_id = 0; spec_id < spec_num_ - 1; spec_id++){
    PeakRow cur_row = matrix_[spec_id];
    PeakRow next_row = matrix_[spec_id+1];
    find_pair_neighbors(cur_row, next_row, search_bin_num, mass_tol);
  }
}

void toppic::PeakMatrix::remove_peak(ExpPeak peak){

  int spec_id = peak.getSpecId();
  int bin_idx = get_index(peak.getPos());
  std::vector<std::vector<ExpPeak>> peaks_row = matrix_[spec_id].getRow();
  std::vector<ExpPeak> bin_peaks = peaks_row[bin_idx];
  for (std::vector<ExpPeak>::iterator bin_peak = bin_peaks.begin(); bin_peak < bin_peaks.end(); bin_peak++) {
    if (bin_peak->getPeakId() == peak.getPeakId()) {
      bin_peaks.erase(bin_peak);
      break;
    }
  }
}

void toppic::PeakMatrix::remove_peak_in_range(int spec_id, double min_pos, double max_pos) {
  int start_bin_idx = get_index(min_pos);
  int end_bin_idx = get_index(max_pos);
  for (int bin_idx = start_bin_idx; bin_idx < end_bin_idx + 1; bin_idx++) {
    std::vector<std::vector<ExpPeak>> peaks_row = matrix_[spec_id].getRow();
    std::vector<ExpPeak> bin_peaks = peaks_row[bin_idx];
    for (std::vector<ExpPeak>::iterator bin_peak = bin_peaks.begin(); bin_peak < bin_peaks.end(); bin_peak++) {
      if (bin_peak->getPos() >= min_pos and bin_peak->getPos() <= max_pos)
        bin_peaks.erase(bin_peak);
    }
  }
}
