//
// Created by abbash on 10/8/2021.
//

#include "feature_extract.h"


namespace toppic {
namespace feature_extract {

double get_mass(double mz, int charge) {
  double proton = 1.00727;
  return ((mz * charge) - (charge * proton));
}

double get_mz(double mass, int charge) {
  double proton = 1.00727;
  return ((mass + (charge * proton)) / charge);
}

std::vector<double> get_base_distribution(DeconvPeakPtr spectrum_peak) {
  std::vector<double> base_distribution;
  std::vector<double> distribution = spectrum_peak->getTheoEnvelope();
  for (size_t i = 0; i < distribution.size(); i++)
    base_distribution.push_back(get_mass(distribution[i], spectrum_peak->getCharge()));
  return base_distribution;
}

std::vector<double> get_charge_distribution(std::vector<double> base_distribution, int charge) {
  std::vector<double> distribution;
  for (size_t i = 0; i < base_distribution.size(); i++)
    distribution.push_back(get_mz(base_distribution[i], charge));
  return distribution;
}

bool chek_output(pair<int, int> feature_boundary, std::vector<double> xic_data) {
  if (std::get<0>(feature_boundary) == -1) {
    LOG_ERROR("********* Null Output!!! ");
    return false;
  }
  double peak_sum = 0.0;
  for (size_t i = std::get<0>(feature_boundary); i <= std::get<1>(feature_boundary); i++)
    peak_sum = peak_sum + xic_data[i];
  if (peak_sum == 0) {
    LOG_ERROR("$$$$$$$$$$$ Not enough peak data!!!");
    return false;
  }
  return true;
}

bool check_feature_length(pair<int, int> refined_feature_boundary, std::vector<double> xic_data) {
  if (std::get<0>(refined_feature_boundary) == -1) {
    LOG_ERROR("Null Output!!! ");
    return false;
  }
  double peak_sum = 0.0;
  for (size_t i = std::get<0>(refined_feature_boundary); i <= std::get<1>(refined_feature_boundary); i++)
    peak_sum = peak_sum + xic_data[i];
  if (peak_sum == 0) {
    LOG_ERROR("Not enough peak data!!!");
    return false;
  }
  if (std::get<1>(refined_feature_boundary) - std::get<0>(refined_feature_boundary) + 1 <= 0) {
    LOG_ERROR("Skipping Feature of length " << 0);
    return false;
  }
  return true;
}

void getSampleFeatures(SampleFeaturePtrVec &sample_features, FracFeaturePtrVec &frac_features) {
  //sample features;
  for (size_t i = 0; i < frac_features.size(); i++) {
    SampleFeaturePtr sample_feature = std::make_shared<SampleFeature>(frac_features[i], frac_features[i]->getId());
    sample_features.push_back(sample_feature);
    frac_features[i]->setSampleFeatureId(frac_features[i]->getId());
    frac_features[i]->setSampleFeatureInte(frac_features[i]->getIntensity());
  }
}

SingleChargeFeaturePtr
return_feature(TopfdParaPtr topfd_para_ptr, mzmlProcessorPtr mzml_processor, SpectralProcessorPtr spectral_processor,
               DeconvPeakPtr spectrum_peak, int peak_idx, int feat_id, pair<int, int> parent_feature_boundary,
               PythonCallerPtr pythonCallerPtr) {
  DeconvPeakPtrVec matched_peaks = spectral_processor->getMatchedPeaks(spectrum_peak, 15E-6, topfd_para_ptr);
  XicGeneratorPtr xic_generator = std::make_shared<XicGenerator>(spectrum_peak, mzml_processor, pythonCallerPtr);
  pair<int, int> feature_boundary = xic_generator->get_feature_boundary(mzml_processor, spectrum_peak,
                                                                        parent_feature_boundary);
  if (!chek_output(feature_boundary, xic_generator->get_xic()))
    return nullptr;

  pair<int, int> refined_feature_boundary = xic_generator->refine_feature_boundary(spectral_processor, feature_boundary,
                                                                                   spectrum_peak, matched_peaks, 3);

  if (!check_feature_length(feature_boundary, xic_generator->get_xic()))
    return nullptr;

  DeconvPeakPtrVec shortlistMatchedPeaks = spectral_processor->shortlistPeaks(matched_peaks,
                                                                              std::get<0>(refined_feature_boundary),
                                                                              std::get<1>(refined_feature_boundary));
  spectral_processor->removePeaks(shortlistMatchedPeaks);
  if (std::get<0>(refined_feature_boundary) == -1)
    xic_generator->remove_peak_data(mzml_processor, refined_feature_boundary);
  else
    xic_generator->remove_peak_data(mzml_processor, parent_feature_boundary);

  SingleChargeFeaturePtr single_feature = std::make_shared<SingleChargeFeature>(spectrum_peak->getCharge(),
                                                                                mzml_processor->get_retention_time(
                                                                                    std::get<0>(
                                                                                        refined_feature_boundary)),
                                                                                mzml_processor->get_retention_time(
                                                                                    std::get<1>(
                                                                                        refined_feature_boundary)),
                                                                                std::get<0>(refined_feature_boundary),
                                                                                std::get<1>(refined_feature_boundary),
                                                                                spectrum_peak->getIntensity(),
                                                                                matched_peaks.size());
  return single_feature;
}

SingleChargeFeaturePtrVec
explore_neighboring_charge_states(TopfdParaPtr topfd_para_ptr, SingleChargeFeaturePtr feature_ptr,
                                  mzmlProcessorPtr mzml_processor, SpectralProcessorPtr spectral_processor,
                                  DeconvPeakPtr spectrum_peak, int peak_idx, int feat_id,
                                  PythonCallerPtr pythonCallerPtr) {
  std::vector<double> base_distribution = get_base_distribution(spectrum_peak);
  SingleChargeFeaturePtrVec neighborFeatureList;
  int miss_count = 0;
  for (int charge = spectrum_peak->getCharge() - 1; charge > -1; charge--) {
    if (miss_count == 2)
      break;
    LOG_ERROR("");
    LOG_ERROR("Processing neighboring charge: " << charge);
    std::vector<double> distribution = get_charge_distribution(base_distribution, charge);
    DeconvPeakPtr temp_peak = std::make_shared<DeconvPeak>(spectrum_peak->getSpId(), spectrum_peak->getId(),
                                                           spectrum_peak->getMonoMass(), spectrum_peak->getIntensity(),
                                                           charge);
    temp_peak->setTheoEnvelope(distribution);
    temp_peak->setTheoEnvelopeInte(spectrum_peak->getTheoEnvelopeInte());
//    LOG_ERROR("GET NEIGHBORING FEATURES " << feature_ptr->getScanBegin() << " " << feature_ptr->getScanEnd());
    SingleChargeFeaturePtr temp_ms1_feature = return_feature(topfd_para_ptr, mzml_processor, spectral_processor,
                                                             temp_peak, peak_idx, feat_id,
                                                             std::make_pair<int, int>(feature_ptr->getScanBegin(),
                                                                                      feature_ptr->getScanEnd()),
                                                             pythonCallerPtr);
    if (temp_ms1_feature == nullptr) {
      miss_count = miss_count + 1;
      continue;
    }
    neighborFeatureList.push_back(temp_ms1_feature);
  }

  miss_count = 0;
  for (int charge = spectrum_peak->getCharge() + 1; topfd_para_ptr->max_charge_; charge++) {
    if (miss_count == 2)
      break;
    LOG_ERROR("");
    LOG_ERROR("Processing neighboring charge: " << charge);
    std::vector<double> distribution = get_charge_distribution(base_distribution, charge);
    DeconvPeakPtr temp_peak = std::make_shared<DeconvPeak>(spectrum_peak->getSpId(), spectrum_peak->getId(),
                                                           spectrum_peak->getMonoMass(), spectrum_peak->getIntensity(),
                                                           charge);
    temp_peak->setTheoEnvelope(distribution);
    temp_peak->setTheoEnvelopeInte(spectrum_peak->getTheoEnvelopeInte());
    SingleChargeFeaturePtr temp_ms1_feature = return_feature(topfd_para_ptr, mzml_processor, spectral_processor,
                                                             temp_peak, peak_idx, feat_id,
                                                             std::make_pair<int, int>(feature_ptr->getScanBegin(),
                                                                                      feature_ptr->getScanEnd()),
                                                             pythonCallerPtr);

    if (temp_ms1_feature == nullptr) {
      miss_count = miss_count + 1;
      continue;
    }
    neighborFeatureList.push_back(temp_ms1_feature);
  }
  return neighborFeatureList;
}

FracFeaturePtr
set_FracFeaturePtr(int feat_id, DeconvPeakPtr spectrum_peak, SingleChargeFeaturePtrVec neighborFeatureList) {
  std::vector<int> scanBegin;
  std::vector<int> scanEnd;
  std::vector<double> TimeBegin;
  std::vector<double> TimeEnd;
  std::vector<int> Charge;

  for (int i = 0; i < neighborFeatureList.size(); i++) {
    scanBegin.push_back(neighborFeatureList[i]->getScanBegin());
    scanEnd.push_back(neighborFeatureList[i]->getScanEnd());
    TimeBegin.push_back(neighborFeatureList[i]->getTimeBegin());
    TimeEnd.push_back(neighborFeatureList[i]->getTimeEnd());
    Charge.push_back(neighborFeatureList[i]->getCharge());
  }

//  FracFeature(int id, int fraction_id,
//  const std::string &file_name,
//  double mono_mass, double inte,
//  int min_ms1_id, int max_ms1_id,
//  double retent_begin, double retent_end,
//  int scan_begin, int scan_end,
//  int min_charge, int max_charge,
//  int env_num, double time_apex);
//
  FracFeaturePtr feature = std::make_shared<FracFeature>(feat_id,
                                                         0, "", spectrum_peak->getMonoMass(),
                                                         spectrum_peak->getIntensity(),
                                                         *min_element(scanBegin.begin(), scanBegin.end()),
                                                         *max_element(scanEnd.begin(), scanEnd.end()),
                                                         *min_element(TimeBegin.begin(), TimeBegin.end()),
                                                         *max_element(TimeEnd.begin(), TimeEnd.end()),
                                                         -1, -1,
                                                         *min_element(Charge.begin(), Charge.end()),
                                                         *max_element(Charge.begin(), Charge.end()),
                                                         neighborFeatureList[0]->getEnvNum(), 0);
  return feature;
}

int process(TopfdParaPtr topfd_para_ptr, std::vector<std::string> spec_file_lst) {
  mzmlProcessorPtr mzml_processor = std::make_shared<mzmlProcessor>(spec_file_lst[0], topfd_para_ptr);

  SpectralProcessorPtr spectral_processor = std::make_shared<SpectralProcessor>(spec_file_lst[1], topfd_para_ptr,
                                                                                mzml_processor);
  LOG_ERROR("Number of peaks: " << spectral_processor->getNumPeaks());

  mzml_processor->filter_peaks_RT(); ///////////////////////////////
  mzml_processor->generate_matrix();

  PythonCallerPtr pythonCallerPtr = std::make_shared<PythonCaller>();

  FracFeaturePtrVec frac_features;
  int feat_id = 0;
  for (int peak_idx = 0; peak_idx < spectral_processor->getNumPeaks(); peak_idx++) {
    LOG_ERROR("\n###################################################### " << peak_idx << " | " << feat_id);
    DeconvPeakPtr spectrum_peak = spectral_processor->get_spec_peak(peak_idx);
    if (spectral_processor->exists(spectrum_peak))
      continue;
    SingleChargeFeaturePtr feature_ptr = return_feature(topfd_para_ptr, mzml_processor, spectral_processor,
                                                        spectrum_peak,
                                                        peak_idx, feat_id, std::make_pair(-1, -1), pythonCallerPtr);
    if (feature_ptr == nullptr)
      continue;
    SingleChargeFeaturePtrVec neighborFeatureList = explore_neighboring_charge_states(topfd_para_ptr, feature_ptr,
                                                                                      mzml_processor,
                                                                                      spectral_processor,
                                                                                      spectrum_peak, peak_idx, feat_id,
                                                                                      pythonCallerPtr);

    neighborFeatureList.push_back(feature_ptr);
    FracFeaturePtr feature = set_FracFeaturePtr(feat_id, spectrum_peak, neighborFeatureList);
    feature->setSingleFeatures(neighborFeatureList);
//    LOG_ERROR("Neighboring Features " << neighborFeatureList.size());
    frac_features.push_back(feature);
    feat_id++;
  }
  LOG_ERROR("");
  LOG_ERROR("Write Output File");
  std::string batmass_file_name = "f_frac.mzrt.csv";
  frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features);

  SampleFeaturePtrVec sample_features;
  getSampleFeatures(sample_features, frac_features);
  std::string sample_feature_file_name = "f_ms1.feature";
  sample_feature_writer::writeFeatures(sample_feature_file_name, sample_features);
  pythonCallerPtr->close_connection();
  return 0;
}

}
}