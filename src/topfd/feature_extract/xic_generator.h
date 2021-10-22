//
// Created by abbash on 9/24/2021.
//

#ifndef TOPPIC_XIC_GENERATOR_H
#define TOPPIC_XIC_GENERATOR_H

#include "ms/spec/deconv_peak.hpp"
#include "common/util/logger.hpp"
#include "topfd/feature_extract/spectral_data_processor.h"
#include "topfd/feature_extract/mzml_data_processor.h"
#include "python_caller.h"

namespace toppic {

class XicGenerator {
 public:
  XicGenerator(DeconvPeakPtr spectrum_peak, mzmlProcessorPtr mzml_processor, PythonCallerPtr pythonCallerPtr);

  pair<int, int> get_feature_boundary(mzmlProcessorPtr mzml_processor, DeconvPeakPtr spectrum_peak,
                                      pair<int, int> parent_feature_boundary);

  pair<int, int> refine_feature_boundary(SpectralProcessorPtr spectral_processor, pair<int, int> feature_boundary,
                                         DeconvPeakPtr spectrum_peak,
                                         DeconvPeakPtrVec matched_peaks, int miss_tol);

  void remove_peak_data(mzmlProcessorPtr mzml_processor, std::pair<int, int> boundary);

  std::vector<double> get_xic() { return data_matrix_; };

 private:
  PythonCallerPtr pythonCallerPtr_;
  PeakPtrVec2D y_exp_data_;
  std::vector<double> data_matrix_;
  std::vector<double> y_fit_;
  double sigma_;

  double _sum_theo_envelope_inte(DeconvPeakPtr spectrum_peak);

  vector<double> _get_normalized_theo_envelope_inte(DeconvPeakPtr spectrum_peak);

  vector<double> _scaleTheoInte(DeconvPeakPtr spectrum_peak, vector<double> exp_inte_distribution);

  bool _evaluate_data_matrix(vector<double> data_matrix, vector<double> normalizedTheoEnvelopeInte, double tol);

  vector<double> _get_normalized_value(vector<double> vec, double normFactor);

  double _pearsonCorrelationCoefficient(vector<double> theoEnvelopeInte, vector<double> expEnvelopeInte);

  pair<int, int> _get_xic_idices(mzmlProcessorPtr mzml_processor, vector<double> y,
                                 DeconvPeakPtr spectrum_peak, std::pair<int, int> parent_feature_boundary);

  bool _check_feature(DeconvPeakPtr spectrum_peak);

  void _set_fitted_gaussian_data_removed(DeconvPeakPtr spectrum_peak);

  void _set_fitted_gaussian_data_removed_explore(DeconvPeakPtr spectrum_peak,
                                                 pair<int, int> parent_feature_boundary);

  std::vector<double> _fit_vogit(std::vector<double> xic, int ref_scan);

};

typedef std::shared_ptr<XicGenerator> XicGeneratorPtr;
}

#endif //TOPPIC_XIC_GENERATOR_H
