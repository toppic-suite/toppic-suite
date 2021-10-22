//
// Created by abbash on 9/24/2021.
//

#include "xic_generator.h"

namespace toppic {

XicGenerator::XicGenerator(DeconvPeakPtr spectrum_peak, mzmlProcessorPtr mzml_processor,
                           PythonCallerPtr pythonCallerPtr) {
  pythonCallerPtr_ = pythonCallerPtr;
  LOG_ERROR("Mass: " << spectrum_peak->getMonoMass() << ", SpID:" << spectrum_peak->getSpId() << ", pId:"
                     << spectrum_peak->getId() << ", C:" << spectrum_peak->getCharge() << ", Inte:"
                     << spectrum_peak->getIntensity());
  y_exp_data_ = mzml_processor->get_data_matrix(spectrum_peak);
  std::vector<double> injection_times = mzml_processor->getInjectionTimes();
  data_matrix_ = std::vector<double>(mzml_processor->get_num_spectra(), 0);
  if (_sum_theo_envelope_inte(spectrum_peak) == 0 or y_exp_data_[0].size() < 2)
    // Here, we check if we have an envelope with positive intenseities
    // we also check if we have data coming from the mzml data
    return;
  std::vector<double> normalizedTheoEnvelopeInte = _get_normalized_theo_envelope_inte(spectrum_peak);
  std::vector<double> TheoEnvelope = spectrum_peak->getTheoEnvelope();
  std::vector<double> TheoEnvelopeInte = spectrum_peak->getTheoEnvelopeInte();

  // compute XIC value
  for (size_t i = 0; i < y_exp_data_[0].size(); i++) {
    std::vector<double> data_matrix;
    std::vector<double> data_matrix_masses;
    for (size_t j = 0; j < y_exp_data_.size(); j++) {
      data_matrix_masses.push_back(y_exp_data_[j][i]->getPosition());
      data_matrix.push_back(y_exp_data_[j][i]->getIntensity());
    }

    if (std::accumulate(data_matrix.begin(), data_matrix.end(), 0) == 0.0)
      continue;

    std::vector<double> scaled_theo_inte = _scaleTheoInte(spectrum_peak, data_matrix);
    bool condition = _evaluate_data_matrix(data_matrix, scaled_theo_inte, 0.5);
    double scalled_theo_inte_sum = 0;
    if (condition)
      scalled_theo_inte_sum = std::accumulate(scaled_theo_inte.begin(), scaled_theo_inte.end(), 0.0);
    //data_matrix_[i] = scalled_theo_inte_sum / injection_times[i]; ////////////////////////// GET INJECTION TIMES
    data_matrix_[i] = scalled_theo_inte_sum;

  }
}

double XicGenerator::_sum_theo_envelope_inte(DeconvPeakPtr spectrum_peak) {
  std::vector<double> theoEnvelopeInte = spectrum_peak->getTheoEnvelopeInte();
  return std::accumulate(theoEnvelopeInte.begin(), theoEnvelopeInte.end(), 0.0);
}

std::vector<double> XicGenerator::_get_normalized_theo_envelope_inte(DeconvPeakPtr spectrum_peak) {
  std::vector<double> theoEnvelopeInte = spectrum_peak->getTheoEnvelopeInte();
  double max_inte = *std::max_element(theoEnvelopeInte.begin(), theoEnvelopeInte.end());
  return _get_normalized_value(theoEnvelopeInte, max_inte);
}

std::vector<double> XicGenerator::_get_normalized_value(std::vector<double> vec, double normFactor) {
  std::vector<double> normalizedVec;
  for (size_t i = 0; i < vec.size(); i++) {
    normalizedVec.push_back(vec[i] / normFactor);
  }
  return normalizedVec;
}

// compute the intensity ratio based on the top three peaks
std::vector<double>
XicGenerator::_scaleTheoInte(DeconvPeakPtr spectrum_peak, std::vector<double> exp_inte_distribution) {
  std::vector<double> theoEnvelopeInte = spectrum_peak->getTheoEnvelopeInte();
  std::vector<double> theoEnvelopePosition = spectrum_peak->getTheoEnvelope();
  double inteRatio;
  double theo_sum = 0;
  double obs_sum = 0;

  int refer_idx = std::max_element(theoEnvelopeInte.begin(), theoEnvelopeInte.end()) - theoEnvelopeInte.begin();
  theo_sum = theo_sum + theoEnvelopeInte[refer_idx];
  obs_sum = obs_sum + exp_inte_distribution[refer_idx];
  if (refer_idx - 1 >= 0) {
    theo_sum = theo_sum + theoEnvelopeInte[refer_idx - 1];
    obs_sum = obs_sum + exp_inte_distribution[refer_idx - 1];
  }
  int size = theoEnvelopeInte.size();
  if (refer_idx + 1 < size) {
    theo_sum = theo_sum + theoEnvelopeInte[refer_idx + 1];
    obs_sum = obs_sum + exp_inte_distribution[refer_idx - 1];
  }
  if (theo_sum == 0)
    inteRatio = 1.0;
  else if (obs_sum == 0)
    inteRatio = 1.0;
  else
    inteRatio = obs_sum / theo_sum;
  std::vector<double> scaled_theo_inte;
  for (size_t i = 0; i < theoEnvelopeInte.size(); i++)
    scaled_theo_inte.push_back(theoEnvelopeInte[i] * inteRatio);
  return scaled_theo_inte;
}

double XicGenerator::_pearsonCorrelationCoefficient(std::vector<double> theoEnvelopeInte,
                                                    std::vector<double> expEnvelopeInte) { //
  double theoEnvelopeInteSum = 0, expEnvelopeInteSum = 0, productSum = 0;
  double theoEnvelopeInteSquareSum = 0, expEnvelopeInteSquareSum = 0;
  int numElem = theoEnvelopeInte.size();
  for (size_t i = 0; i < theoEnvelopeInte.size(); i++) { // compute summations
    theoEnvelopeInteSum = theoEnvelopeInteSum + theoEnvelopeInte[i];
    expEnvelopeInteSum = expEnvelopeInteSum + expEnvelopeInte[i];
    productSum = productSum + (theoEnvelopeInte[i] * expEnvelopeInte[i]);
    theoEnvelopeInteSquareSum = theoEnvelopeInteSquareSum + (theoEnvelopeInte[i] * theoEnvelopeInte[i]);
    expEnvelopeInteSquareSum = expEnvelopeInteSquareSum + (expEnvelopeInte[i] * expEnvelopeInte[i]);
  }
  return (numElem * productSum - theoEnvelopeInteSum * expEnvelopeInteSum) / sqrt(
      (numElem * theoEnvelopeInteSquareSum - theoEnvelopeInteSum * theoEnvelopeInteSum) *
      (numElem * expEnvelopeInteSquareSum - expEnvelopeInteSum * expEnvelopeInteSum));
}

bool
XicGenerator::_evaluate_data_matrix(std::vector<double> data_matrix, std::vector<double> normalizedTheoEnvelopeInte,
                                    double tol) {
  double max_exp_inte = *std::max_element(data_matrix.begin(), data_matrix.end());
  std::vector<double> normalized_data_matrix = _get_normalized_value(data_matrix, max_exp_inte);
  double corr = _pearsonCorrelationCoefficient(normalized_data_matrix, normalizedTheoEnvelopeInte);
  if (corr >= tol)
    return true;
  else
    return false;
}

std::pair<int, int>
XicGenerator::_get_xic_idices(mzmlProcessorPtr mzml_processor, std::vector<double> y, DeconvPeakPtr spectrum_peak,
                              std::pair<int, int> parent_feature_boundary) {
  std::vector<double> theoEnvelopeInte = spectrum_peak->getTheoEnvelopeInte();
  double mean_base_inte = mzml_processor->get_mean_base_inte();
  double max_inte = *std::max_element(theoEnvelopeInte.begin(), theoEnvelopeInte.end());
  int center = spectrum_peak->getSpId();
  if (std::get<0>(parent_feature_boundary) > -1)
    center = center - std::get<0>(parent_feature_boundary);

  double limit_inte = (0.1 * mean_base_inte) / max_inte;
  LOG_ERROR("Max Inte: " << max_inte << "   " << mean_base_inte << "  " << limit_inte);
  int start_idx = 0;
  for (int idx = center - 1; idx > -1; idx--) {
    double point = y[idx];
//    double limit_inte = (mzml_processor->get_sp_base_inte(idx) / max_inte);
    if (point <= limit_inte) {
      start_idx = idx;
      break;
    }
  }
  int end_idx = y.size() - 1;
  for (int idx = center + 1; idx < y.size(); idx++) {
    double point = y[idx];
//    double limit_inte = (mzml_processor->get_sp_base_inte(idx) / max_inte);
    if (point <= limit_inte) {
      end_idx = idx;
      break;
    }
  }

  if (std::get<0>(parent_feature_boundary) > -1) {
    start_idx = start_idx + std::get<0>(parent_feature_boundary);
    end_idx = end_idx + std::get<0>(parent_feature_boundary);
  }
  LOG_ERROR("Boundary: (" << start_idx << ", " << end_idx << ")");
  return std::make_pair(start_idx, end_idx);
}

std::pair<int, int> XicGenerator::get_feature_boundary(mzmlProcessorPtr mzml_processor, DeconvPeakPtr spectrum_peak,
                                                       std::pair<int, int> parent_feature_boundary) {
  if (!_check_feature(spectrum_peak)) {
//    LOG_ERROR("Check Failed");
    return std::make_pair(-1, -1);
  }

  if (std::get<0>(parent_feature_boundary) == -1)
    _set_fitted_gaussian_data_removed(spectrum_peak);
  else
    _set_fitted_gaussian_data_removed_explore(spectrum_peak, parent_feature_boundary);
  return _get_xic_idices(mzml_processor, y_fit_, spectrum_peak, parent_feature_boundary);
}

bool XicGenerator::_check_feature(DeconvPeakPtr spectrum_peak) {
  if (std::accumulate(data_matrix_.begin(), data_matrix_.end(), 0) == 0.0) {
    return false;
  }

//  double neighbors = 0;
//  for (size_t idx = spectrum_peak->getSpId() - 2; idx <= spectrum_peak->getSpId() + 2; idx++)
//    neighbors = neighbors + data_matrix_[idx];
//  if (neighbors <= 1) {
//    LOG_ERROR("neighbor caused -- ");
//    return false;
//  }
  return true;
}

void XicGenerator::_set_fitted_gaussian_data_removed(DeconvPeakPtr spectrum_peak) {
  double norm_factor = *std::max_element(data_matrix_.begin(), data_matrix_.end());
  vector<double> noramlized_xic;
  for (size_t idx = 0; idx < data_matrix_.size(); idx++) {
    noramlized_xic.push_back(data_matrix_[idx] / norm_factor);
  }
  y_fit_ = pythonCallerPtr_->process(noramlized_xic, spectrum_peak->getSpId());
}

void XicGenerator::_set_fitted_gaussian_data_removed_explore(DeconvPeakPtr spectrum_peak,
                                                             std::pair<int, int> parent_feature_boundary) {

  int start = std::get<0>(parent_feature_boundary);
  int end = std::get<1>(parent_feature_boundary);
  while (end - start + 1 <= 4) {
    int size = data_matrix_.size();
    if (end < size and start > 0) {
      start = start - 1;
      end = end + 1;
    } else if (start == 0) {
      start = start;
      end = end + 2;
    } else if (end == size) {
      start = start - 2;
      end = end;
    }
  }
  parent_feature_boundary = std::make_pair(start, end);
  double norm_factor = *std::max_element(data_matrix_.begin(), data_matrix_.end());
  vector<double> noramlized_xic;
  for (int idx = start; idx <= end; idx++)
    noramlized_xic.push_back(data_matrix_[idx] / norm_factor);
  y_fit_ = pythonCallerPtr_->process(noramlized_xic, spectrum_peak->getSpId() - start);
}

std::pair<int, int>
XicGenerator::refine_feature_boundary(SpectralProcessorPtr spectral_processor, std::pair<int, int> feature_boundary,
                                      DeconvPeakPtr spectrum_peak, DeconvPeakPtrVec matched_peaks, int miss_tol) {
  int ref_sp_id = spectrum_peak->getSpId();
  DeconvPeakPtrVec shortlistMatchedPeaks = spectral_processor->shortlistPeaks(matched_peaks, 0, data_matrix_.size());
  std::vector<int> spectra_id;
  for (size_t idx = 0; idx < shortlistMatchedPeaks.size(); idx++)
    spectra_id.push_back(shortlistMatchedPeaks[idx]->getSpId());
  int refined_start = ref_sp_id;
  int miss_num = 0;
  while (refined_start > 0) {
    if (data_matrix_[refined_start] != 0) {
      if (std::count(spectra_id.begin(), spectra_id.end(), refined_start))
        miss_num = 0;
    } else {
      miss_num = miss_num + 1;
    }
    if (miss_num > miss_tol) {
      refined_start = refined_start + miss_tol;
      break;
    }
    refined_start = refined_start - 1;
  }

  int refined_end = ref_sp_id;
  miss_num = 0;
  while (refined_end < data_matrix_.size() - 1) {
    if (data_matrix_[refined_end] != 0) {
      if (std::count(spectra_id.begin(), spectra_id.end(), refined_end))
        miss_num = 0;
    } else {
      miss_num = miss_num + 1;
    }
    if (miss_num > miss_tol) {
      refined_end = refined_end - miss_tol;
      break;
    }
    refined_end = refined_end + 1;
  }
  LOG_ERROR("Refined Boundary: (" << refined_start << ", " << refined_end << ")");
  return std::make_pair<>(refined_start, refined_end);
}

void XicGenerator::remove_peak_data(mzmlProcessorPtr mzml_processor, std::pair<int, int> boundary) {
  for (int rt_idx = std::get<0>(boundary); rt_idx <= std::get<1>(boundary); rt_idx++) {
    if (rt_idx < 0 or rt_idx > mzml_processor->get_num_spectra())
      continue;
    for (size_t elem_idx = 0; elem_idx < y_exp_data_.size(); elem_idx++) {
      PeakPtr elem = y_exp_data_[elem_idx][rt_idx];
      if (elem->getIntensity() == 0 and elem->getPosition() == 0)
        continue;
      PeakPtrVec peaks = mzml_processor->get_spectrum(rt_idx);
      for (size_t j = 0; j < peaks.size(); j++) {
        PeakPtr peak = peaks[j];
        if ((elem->getPosition() - peak->getPosition()) < 0.005) {
          peaks[j]->setPosition(0);
          peaks[j]->setIntensity(0);

        }
      }
    }
  }
}

} //namespace