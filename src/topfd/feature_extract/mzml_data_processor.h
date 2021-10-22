//
// Created by abbash on 9/23/2021.
//

#ifndef TOPPIC_MZML_DATA_PROCESSOR_H
#define TOPPIC_MZML_DATA_PROCESSOR_H

#include "topfd/common/topfd_para.hpp"
#include "ms/spec/deconv_peak.hpp"
#include "topfd/msreader/raw_ms_reader.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/spec/peak.hpp"
#include "common/util/logger.hpp"

namespace toppic {

class mzmlProcessor {
 public:
  mzmlProcessor(std::string spec_file_name, TopfdParaPtr topfd_para_ptr);

  int get_index(double mz);

  void generate_matrix();

  PeakPtrVec2D get_data_matrix(DeconvPeakPtr spectrum_peak);

  PeakPtrVec get_exp_data(double peak, double tol);

  const vector<double> &getInjectionTimes() const {
    return injection_times_;
  }

  void setInjectionTimes(const vector<double> &injectionTimes) {
    injection_times_ = injectionTimes;
  }

  int get_num_spectra() { return peakData_.size(); }

  PeakPtrVec get_spectrum(int id) { return peakData_[id]; }

  double get_sp_base_inte(int idx) { return base_intensity_[idx]; }

  double get_retention_time(int idx) { return retention_times_[idx]; }

  void filter_peaks_RT();

  double get_mean_base_inte() { return mean_base_inte_; }

 private:
  double round_(double i); // round to 3 digits for indexing

  void _get_exp_data(double peak, double tol);

  void update_peak_location();

  PeakPtrVec2D peakData_;
  std::vector<double> base_intensity_;
  std::vector<double> retention_times_;
  std::vector<double> injection_times_;
  double min_mz_ = std::numeric_limits<double>::max();
  double max_mz_ = std::numeric_limits<double>::min();
  double bin_size_ = 0.1;
  double mean_base_inte_ = 0;
  std::vector<std::vector<std::pair<int, int>>> peak_locations_;


};

typedef std::shared_ptr<mzmlProcessor> mzmlProcessorPtr;
}

#endif //TOPPIC_MZML_DATA_PROCESSOR_H
