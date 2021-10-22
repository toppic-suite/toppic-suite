//
// Created by abbash on 9/17/2021.
//

#include "spectral_data_processor.h"

namespace toppic {

SpectralProcessor::SpectralProcessor(std::string ms1_file_name, TopfdParaPtr topfd_para_ptr,
                                     mzmlProcessorPtr mzml_processor) {
  LOG_ERROR(" Spectral Processing ");
  EnvBase::initBase(topfd_para_ptr->resource_dir_);
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>(topfd_para_ptr);
  SimpleMsAlignReader::readMsOneSpectra(ms1_file_name, spec_list_);
  std::vector<double> intes;
  for (size_t i = 0; i < spec_list_.size(); i++) {
    DeconvPeakPtrVec peaks = spec_list_[i]->getPeakPtrVec();
    PeakPtrVec peak_list = mzml_processor->get_spectrum(i);
    for (size_t base = 0; base < peaks.size(); base++) {
      get_envelope(peaks, base, peak_list, env_para_ptr);
    }
    spectrum_peaks_.insert(std::end(spectrum_peaks_), std::begin(peaks), std::end(peaks));
    spectrum_peaks_mass_sorted_.insert(std::end(spectrum_peaks_mass_sorted_), std::begin(peaks),
                                       std::end(peaks));
  }
  std::sort(spectrum_peaks_.begin(), spectrum_peaks_.end(), Peak::cmpInteDec);
  std::sort(spectrum_peaks_mass_sorted_.begin(), spectrum_peaks_mass_sorted_.end(), Peak::cmpMassDec);
}

void
SpectralProcessor::get_envelope(DeconvPeakPtrVec peaks, size_t base, PeakPtrVec peak_list, EnvParaPtr env_para_ptr) {
  double base_mass = peaks[base]->getPosition();
  int charge = peaks[base]->getCharge();
  double base_mz = peaks[base]->getMonoMz();

  EnvelopePtr ref_env = EnvBase::getStaticEnvByBaseMass(base_mass);
  EnvelopePtr theo_env = ref_env->distrToTheoMono(base_mz, charge);

  double ratio = env_detect::calcInteRatio(peak_list, theo_env, env_para_ptr->getMzTolerance());
  theo_env->changeIntensity(ratio);

  int mass_group = env_para_ptr->getMassGroup(base_mass);
  double percentage = env_para_ptr->getPercentBound(mass_group);
  theo_env = theo_env->getSubEnv(percentage, env_para_ptr->min_refer_inte_,
                                 env_para_ptr->max_back_peak_num_, env_para_ptr->max_forw_peak_num_);
  RealEnvPtr real_env = std::make_shared<RealEnv>(peak_list, theo_env, env_para_ptr->getMzTolerance(),
                                                  env_para_ptr->min_refer_inte_);
  MatchEnvPtr match_env = std::make_shared<MatchEnv>(mass_group, theo_env, real_env);

  theo_env = match_env->getTheoEnvPtr();
  std::vector<double> dist;
  std::vector<double> dist_inte;
  for (int p_idx = 0; p_idx < theo_env->getPeakNum(); p_idx++) {
    dist.push_back(theo_env->getMz(p_idx));
    dist_inte.push_back(theo_env->getIntensity(p_idx));
  }
  peaks[base]->setTheoEnvelope(dist);
  peaks[base]->setTheoEnvelopeInte(dist_inte);
}

int SpectralProcessor::_binary_search(double mass) {
  int low = 0;
  int mid = 0;
  int high = spectrum_peaks_mass_sorted_.size() - 1;
  while (low <= high) {
    mid = (high + low) / 2;
    if (spectrum_peaks_mass_sorted_[mid]->getMonoMass() < mass)
      low = mid + 1;
    else if (spectrum_peaks_mass_sorted_[mid]->getMonoMass() > mass)
      high = mid - 1;
  }
  return low;
}

std::vector<double> SpectralProcessor::_getExtendMasses(double mass) {
  double IM = 1.00235;
  std::vector<double> extend_offsets_{0, -IM, IM};
  std::vector<double> result;
  for (size_t i = 0; i < extend_offsets_.size(); i++) {
    double new_mass = mass + extend_offsets_[i];
    result.push_back(new_mass);
  }
  return result;
}

bool SpectralProcessor::_check_overlap(DeconvPeakPtr spectrum_peak, DeconvPeakPtr peak, TopfdParaPtr topfd_para_ptr) {
  std::tuple<double, double> overlap;
  LOG_ERROR("INCLUDE IMPLEMENTATION");
  return true;
}

DeconvPeakPtrVec
SpectralProcessor::getMatchedPeaks(DeconvPeakPtr spectrum_peak, double tol, TopfdParaPtr topfd_para_ptr) {
  double prec_mass = spectrum_peak->getMonoMass();
  int charge = spectrum_peak->getCharge();
  double error_tole = prec_mass * tol;
  std::vector<double> ext_masses = _getExtendMasses(prec_mass);

  double min_mass = ext_masses[0] - error_tole;
  double max_mass = ext_masses[0] + error_tole;
  if (ext_masses.size() > 1) {
    min_mass = ext_masses[1] - error_tole;
    max_mass = ext_masses[ext_masses.size() - 1] + error_tole;
  }
  int start_idx = _binary_search(min_mass);
  int end_idx = _binary_search(max_mass);
  DeconvPeakPtrVec tmp_matched_peaks(spectrum_peaks_mass_sorted_.begin() + start_idx,
                                     spectrum_peaks_mass_sorted_.begin() + end_idx);

  DeconvPeakPtrVec matched_peaks;
  for (size_t peak_idx = 0; peak_idx < tmp_matched_peaks.size(); peak_idx++) {
    DeconvPeakPtr peak = tmp_matched_peaks[peak_idx];
//  if (peak->getCharge() == charge and _check_overlap(spectrum_peak, peak, topfd_para_ptr)) {
    if (peak->getCharge() == charge) {
      DeconvPeakPtrVec peaks = spec_list_[peak->getSpId()]->getPeakPtrVec();
      if (peaks[peak->getId()]->getId() > -1)
        matched_peaks.push_back(peak);
    }
  }
  std::sort(matched_peaks.begin(), matched_peaks.end(), DeconvPeak::cmpPosInc);
  return matched_peaks;
}

DeconvPeakPtrVec
SpectralProcessor::shortlistPeaks(DeconvPeakPtrVec matched_peaks, int rt_begin_sp_id, int rt_end_sp_id) {
  DeconvPeakPtrVec shortlistMatchedPeaks;
  for (size_t peak_idx = 0; peak_idx < matched_peaks.size(); peak_idx++) {
    DeconvPeakPtr peak = matched_peaks[peak_idx];
    if (peak->getSpId() >= rt_begin_sp_id and peak->getSpId() <= rt_end_sp_id) {
      shortlistMatchedPeaks.push_back(peak);
    }
  }
  return shortlistMatchedPeaks;
}

void SpectralProcessor::removePeaks(DeconvPeakPtrVec shortlistMatchedPeaks) {
  for (size_t i = 0; i < shortlistMatchedPeaks.size(); i++) {
    DeconvPeakPtr peak = shortlistMatchedPeaks[i];
    int sp_id = peak->getSpId();
    int peak_id = peak->getId();
    spec_list_[sp_id]->getPeakPtr(peak_id)->setUsedStatus(true);
  }
}

bool SpectralProcessor::exists(DeconvPeakPtr peak) {
  int sp_id = peak->getSpId();
  int peak_id = peak->getId();
  DeconvPeakPtr remain_peak = spec_list_[sp_id]->getPeakPtr(peak_id);
  return remain_peak->getUsedStatus();
}
} // namespace
