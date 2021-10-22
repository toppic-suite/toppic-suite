//
// Created by abbash on 9/23/2021.
//

#include "mzml_data_processor.h"

namespace toppic {

mzmlProcessor::mzmlProcessor(std::string spec_file_name, TopfdParaPtr topfd_para_ptr) {
  LOG_ERROR("Read mzML data");
  PwMsReaderPtr reader_ptr_ = std::make_shared<PwMsReader>(spec_file_name);
  while (true) {
    reader_ptr_->readNext();
    MsHeaderPtr header_ptr = reader_ptr_->getHeaderPtr();
    if (header_ptr == nullptr) {
      break;
    }
    if (header_ptr->getMsLevel() == 1) {
      // compute baseline intensity (noise intensity level)
      PeakPtrVec peak_list = reader_ptr_->getPeakList();
      std::vector<double> intes;
      for (size_t i = 0; i < peak_list.size(); i++) {
        intes.push_back(peak_list[i]->getIntensity());
      }
      double base_inte = baseline_util::getBaseLine(intes);
      // filter data based on noise intensity level
      double min_refer_inte_ = base_inte * topfd_para_ptr->ms_one_sn_ratio_;
      PeakPtrVec filtered_peak_list;
      for (size_t i = 0; i < peak_list.size(); i++) {
        PeakPtr peak = peak_list[i];
        if (peak->getIntensity() >= min_refer_inte_)
          filtered_peak_list.push_back(peak);
      }

      std::sort(filtered_peak_list.begin(), filtered_peak_list.end(), Peak::cmpMassDec);
      // assign values
      if (min_mz_ > filtered_peak_list[0]->getPosition())
        min_mz_ = filtered_peak_list[0]->getPosition();
      if (max_mz_ < filtered_peak_list[filtered_peak_list.size() - 1]->getPosition())
        max_mz_ = filtered_peak_list[filtered_peak_list.size() - 1]->getPosition();
      peakData_.push_back(filtered_peak_list);
      base_intensity_.push_back(base_inte);
      retention_times_.push_back(header_ptr->getRetentionTime());
//      injection_times_.push_back(header_ptr->getInjectionTime());
    }
  }
  double mean_base_inte = 0;
  for (int i = 0; i < get_num_spectra(); i++)
    mean_base_inte = mean_base_inte + base_intensity_[i];
  mean_base_inte_ = mean_base_inte / get_num_spectra();
}

void mzmlProcessor::filter_peaks_RT() {
  LOG_ERROR("Filter MZ data");
  for (size_t idx = 0; idx < peakData_.size() - 1; idx++) {
    for (size_t peakIdx = 0; peakIdx < peakData_[idx].size(); peakIdx++) {
      PeakPtr peak = peakData_[idx][peakIdx];
      int shortlistedPeakIdx = -1;
      for (size_t neighborPeakIdx = 0; neighborPeakIdx < peakData_[idx + 1].size(); neighborPeakIdx++) {
        PeakPtr neighborPeak = peakData_[idx + 1][neighborPeakIdx];
        if (neighborPeak->getNeighbor() == false) {
          double peak_diff = std::numeric_limits<double>::max();
          double mass_diff = abs(neighborPeak->getPosition() - peak->getPosition());
          if (mass_diff < 0.01 and mass_diff < peak_diff) {
            peak_diff = mass_diff;
            shortlistedPeakIdx = neighborPeakIdx;
          }
          if (neighborPeak->getPosition() - peak->getPosition() > 0.1) {
            if (shortlistedPeakIdx > -1) {
              peak->setNeighbor(true);
              peakData_[idx + 1][shortlistedPeakIdx]->setNeighbor(true);
              break;
            }
          }
        }
      }
    }
  }

  PeakPtrVec2D shortlisted_peaks;
  for (size_t idx = 0; idx < peakData_.size(); idx++) {
    PeakPtrVec temp_peaks;
    for (size_t peakIdx = 0; peakIdx < peakData_[idx].size(); peakIdx++) {
      if (peakData_[idx][peakIdx]->getNeighbor() == true)
        temp_peaks.push_back(peakData_[idx][peakIdx]);
    }
    shortlisted_peaks.push_back(temp_peaks);
  }
  peakData_ = shortlisted_peaks;

}

int mzmlProcessor::get_index(double mz) {
  double mz_diff = mz - min_mz_;
  return int(mz_diff / bin_size_);
}

void mzmlProcessor::generate_matrix() {
  LOG_ERROR("Generating matrix")
  int num_bins = int((max_mz_ - min_mz_) / bin_size_) + 1;
  peak_locations_ = std::vector<std::vector<std::pair<int, int>>>(num_bins,
      std::vector<std::pair<int, int>>(peakData_.size(),
      std::make_pair(-1,
                     -1)));
  for (size_t i = 0; i < peakData_.size(); i++) {
    PeakPtrVec peak_list = peakData_[i];
    int j_start = 0;
    for (size_t bin_idx = 0; bin_idx < num_bins; bin_idx++) {
      double bin_min_mz = round_(min_mz_ + (bin_idx * bin_size_));
      double bin_max_mz = round_(min_mz_ + ((bin_idx + 1) * bin_size_));
      int start = j_start;
      int end = -1;
      for (size_t j = j_start; j < peak_list.size(); j++) {
        double mz = peak_list[j]->getPosition();
        if (mz < bin_min_mz)
          continue;
        if (mz >= bin_max_mz) {
          j_start = j;
          break;
        }
        end = j;
      }
      if (end != -1)
        peak_locations_[bin_idx][i] = std::make_pair(start, end);
    }
  }
}

PeakPtrVec mzmlProcessor::get_exp_data(double peak, double tol) {
  int start_matrix_idx = get_index(round_(peak - tol));
  int end_matrix_idx = get_index(round_(peak + tol));
  std::vector<std::vector<std::pair<int, int>>> temp_peak_location;
  for (size_t idx = start_matrix_idx; idx < end_matrix_idx + 1; idx++)
    temp_peak_location.push_back(peak_locations_[idx]);

  PeakPtrVec out_data;
  for (size_t sp_idx = 0; sp_idx < temp_peak_location[0].size(); sp_idx++) {
    PeakPtrVec peak_list = peakData_[sp_idx];

    std::vector<std::pair<int, int>> location;
    for (size_t i = 0; i < temp_peak_location.size(); i++)
      location.push_back(temp_peak_location[i][sp_idx]);

    PeakPtrVec value;
    for (size_t i = 0; i < location.size(); i++) {
      std::pair<int, int> elem = location[i];
      if (std::get<0>(elem) == -1)
        continue;
      for (size_t k = std::get<0>(elem); k < std::get<1>(elem) + 1; k++)
        value.push_back(peak_list[k]);
    }
    PeakPtr out_val = std::make_shared<Peak>(0, 0);
    if (value.size() == 1) {
      if (std::abs(value[0]->getPosition() - peak) <= tol) {
        out_val->setIntensity(value[0]->getIntensity());
        out_val->setPosition(value[0]->getPosition());
      }
    } else {
      double max_inte = std::numeric_limits<double>::min();
      for (size_t i = 0; i < value.size(); i++) {
        PeakPtr val = value[i];
        if (std::abs(val->getPosition() - peak) <= tol) {
          if (val->getIntensity() > max_inte) {
            max_inte = val->getIntensity();
            out_val->setIntensity(value[i]->getIntensity());
            out_val->setPosition(value[i]->getPosition());
          }
        }
      }
    }
    out_data.push_back(out_val);
  }
  return out_data;
}

PeakPtrVec2D mzmlProcessor::get_data_matrix(DeconvPeakPtr spectrum_peak) {
  std::vector<double> distribution = spectrum_peak->getTheoEnvelope();
  std::vector<double> distribution_inte = spectrum_peak->getTheoEnvelopeInte();

  PeakPtrVec2D data_matrix;
  std::vector<double> new_dist;
  std::vector<double> new_dist_intes;
  for (size_t distribution_peak_idx = 0; distribution_peak_idx < distribution.size(); distribution_peak_idx++) {
    double theo_peak_mass = distribution[distribution_peak_idx];
    if (theo_peak_mass < min_mz_ or theo_peak_mass > max_mz_)
      continue;
    PeakPtrVec exp_dist = get_exp_data(theo_peak_mass, 0.005);
    data_matrix.push_back(exp_dist);

    new_dist.push_back(theo_peak_mass);
    new_dist_intes.push_back(distribution_inte[distribution_peak_idx]);
  }
  spectrum_peak->setTheoEnvelope(new_dist);
  spectrum_peak->setTheoEnvelopeInte(new_dist_intes);
  return data_matrix;
}

double mzmlProcessor::round_(double i) {
  int val = i * 1000.0;
  double rounded_val = val / 1000.0 + 0.001;
  return rounded_val;
}

}