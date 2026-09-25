// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
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

#include "topfd/deconv/flashdeconv_ms2_process.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "ms/env/avg_mass_util.hpp"
#include "ms/env/env.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/env_para.hpp"
#include "ms/env/exp_env.hpp"
#include "ms/env/match_env.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/feature/spec_feature_reader.hpp"
#include "ms/mzml/mzml_ms_group_reader.hpp"
#include "ms/mzml/mzml_ms_sql_writer.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/spec/msalign_thread_merge.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/peak_util.hpp"

namespace toppic {

// The FLASHDeconv MS/MS deconvolution, ported from OpenMS (Jeong et al.,
// FLASHDeconv; src/topp/FLASHDeconv.cpp with ANALYSIS/TOPDOWN/
// SpectralDeconvolution, PeakGroup, PeakGroupScoring and FLASHHelperClasses).
// The steps for one spectrum, each a function below with the OpenMS name it
// comes from in its comment:
//   1. keep the strongest peaks (FLASHDeconvAlgorithm::filterLowPeaks_),
//   2. put every peak on a logarithmic m/z axis (updateLogMzPeaks_),
//   3. bin the log m/z axis with a bin of tolerance / 2.5 (binLogMzPeaks_),
//   4. read the charge-mass matrix through the "universal pattern"
//      log(mass) = log(m/z - proton) + log(z): a mass bin is a candidate when
//      consecutive charge states of similar intensity, or an isotope peak at
//      a low charge, support it and no harmonic peak sits between
//      (updateCandidateMassBins_), and every m/z bin keeps its two strongest
//      masses (filterMassBins_),
//   5. collect the isotope peaks of every candidate mass at every charge into
//      a peak group (getCandidatePeakGroups_),
//   6. refine every peak group: recruit the raw peaks again for the current
//      monoisotopic mass, trim the charge range, score the charge
//      distribution, the averagine isotope cosine (with a monoisotopic
//      offset search), the per-charge cosine, the signal-to-noise ratio and
//      the logistic QScore (PeakGroup::updateQscore, PeakGroupScoring),
//   7. remove charge-error, harmonic and duplicate masses and keep the peak
//     groups above the QScore and SNR thresholds (scoreAndFilterPeakGroups_).
// The theoretical isotope distributions are TopPIC's averagine table
// (EnvBase, res/base_data/theo_patt.txt) prepared as FLASHDeconv prepares its
// own (trimmed to 99.99% of the power and L2-normalised, see Averagine).
namespace flashdeconv_ms2_process {

// ----------------------------------------------------------------------------
// Constants of FLASHDeconv (SpectralDeconvolution.cpp/.h, PeakGroup.h,
// Constants.h)
// ----------------------------------------------------------------------------

// Constants::PROTON_MASS_U; positive ionization mode
constexpr double kProtonMass = 1.0072764667710;
// Constants::ISOTOPE_MASSDIFF_55K_U: the isotope spacing FLASHDeconv uses
constexpr double kIsoDaDistance = 1.002371;
// -SD:tol for MS2 (ppm); FLASHDeconv's own tolerance, not TopFD's -e
constexpr double kMs2TolerancePpm = 10.0;
// tol_div_factor: the log m/z bin is tolerance / 2.5
constexpr double kTolDivFactor = 2.5;
// charges up to low_charge_ are supported by isotope peaks, higher charges by
// charge continuity
constexpr int kLowCharge = 10;
// min_support_peak_count_: charge states needed to make a mass a candidate
constexpr int kMinSupportPeakCount = 2;
// min_iso_size: isotopes needed for a cosine
constexpr int kMinIsoSize = 2;
// min_qscore_, min_cos[MS2], min_snr[MS2]
constexpr double kMinQscore = 0.2;
constexpr double kMinIsotopeCosine = 0.85;
constexpr double kMinSnr = 0.25;
// PeakGroup::updateQscore: a charge distribution score below this is 0
constexpr double kMinChargeScore = 0.7;
// max_mass_dalton_tolerance_: the isotope peak tolerance is at most 0.16 Da
constexpr double kMaxMassDaltonTolerance = 0.16;
// DeconvolvedSpectrum::getCurrentMinMass: 50 Da for every MSn, n > 1
constexpr double kMinMass = 50.0;
// PeakGroup::min_negative_isotope_index_
constexpr int kMinNegativeIsotopeIndex = -1;
// scoreAndFilterPeakGroups_: iterations of the monoisotopic mass refinement
constexpr int kRefineIterationNum = 10;
// FLASHDeconvAlgorithm::filterLowPeaks_ (centroided spectra)
constexpr std::size_t kMaxPeakNum = 30000;
constexpr double kMinPeakInteRatio = 1e-3;
// harmonic_charges_
const std::vector<int> kHarmonicCharges = {2, 3, 5, 7, 11};
// PeakGroupScoring::weight_ (the last one is the intercept)
const std::vector<double> kQscoreWeights = {-21.0476, 1.5045, -0.1303,
                                            0.183,    0.1834, 17.804};
// PrecalculatedAveragine: mass range and step of the table
constexpr double kAveragineMinMass = 50.0;
constexpr double kAveragineMassInterval = 25.0;

// ----------------------------------------------------------------------------
// Averagine (FLASHHelperClasses::PrecalculatedAveragine)
// ----------------------------------------------------------------------------

// The theoretical isotope distributions of masses 50 Da to max_mass in 25 Da
// steps, taken from EnvBase and prepared as PrecalculatedAveragine prepares
// the CoarseIsotopePatternGenerator distributions: the weak end isotopes are
// dropped until 0.01% of the power is gone (at least 2 isotopes stay), the
// intensities are divided by the square root of the remaining power (so that
// the cosine of getCosine needs no norm of the theoretical side), and the
// apex index, the isotope counts left and right of the apex (a running
// maximum over the masses, at least 2), the average-mass offset and the SNR
// multiplication factor (the squared sum of the normalised intensities) are
// recorded.
class Averagine {
 public:
  explicit Averagine(double max_mass);

  const std::vector<double>& getInte(double mass) const {
    return entries_[massToIndex(mass)].inte;
  }

  int getApexIndex(double mass) const {
    return entries_[massToIndex(mass)].apex_index;
  }

  int getLeftCountFromApex(double mass) const {
    return entries_[massToIndex(mass)].left_count;
  }

  int getRightCountFromApex(double mass) const {
    return entries_[massToIndex(mass)].right_count;
  }

  int getLastIndex(double mass) const {
    const Entry& e = entries_[massToIndex(mass)];
    return e.apex_index + e.right_count;
  }

  double getAverageMassDelta(double mass) const {
    return entries_[massToIndex(mass)].avg_mass_delta;
  }

  double getSnrMulFactor(double mass) const {
    return entries_[massToIndex(mass)].snr_mul_factor;
  }

 private:
  struct Entry {
    std::vector<double> inte;
    int apex_index = 0;
    int left_count = 0;
    int right_count = 0;
    double avg_mass_delta = 0;
    double snr_mul_factor = 0;
  };

  std::vector<Entry> entries_;

  std::size_t massToIndex(double mass) const {
    double idx = std::round(std::max(0.0, mass - kAveragineMinMass) /
                            kAveragineMassInterval);
    return std::min(static_cast<std::size_t>(idx), entries_.size() - 1);
  }
};

Averagine::Averagine(double max_mass) {
  const double min_pwr = 0.9999;
  const std::size_t min_iso_length = 2;
  const int min_left_right_count = 2;
  int max_left_count = 0;
  int max_right_count = 0;
  for (double mass = kAveragineMinMass; mass <= max_mass;
       mass += kAveragineMassInterval) {
    EnvPtr env_ptr = EnvBase::getEnvByMonoMass(mass);
    if (env_ptr == nullptr) {
      break;
    }
    // the distribution: mass offsets from the monoisotopic peak and
    // intensities
    int num = env_ptr->getPeakNum();
    std::vector<double> mz(num);
    std::vector<double> inte(num);
    double total_pwr = 0;
    int apex = 0;
    double apex_inte = 0;
    for (int k = 0; k < num; k++) {
      mz[k] = env_ptr->getMz(k);
      inte[k] = env_ptr->getInte(k);
      total_pwr += inte[k] * inte[k];
      if (inte[k] > apex_inte) {
        apex_inte = inte[k];
        apex = k;
      }
    }
    // trim the weaker end until 0.01% of the power is removed
    int left = 0;
    int right = num - 1;
    int trim_count = 0;
    double pwr = 0;
    while (static_cast<std::size_t>(num - trim_count) > min_iso_length &&
           left < right) {
      double lint = inte[left];
      double rint = inte[right];
      bool trim_left = true;
      if (lint < rint) {
        pwr += lint * lint;
      } else {
        pwr += rint * rint;
        trim_left = false;
      }
      if (total_pwr - pwr < total_pwr * min_pwr) {
        break;
      }
      trim_count++;
      if (trim_left) {
        inte[left] = 0;
        left++;
      } else {
        inte[right] = 0;
        right--;
      }
    }
    total_pwr -= pwr;
    Entry entry;
    entry.apex_index = apex;
    int left_count = std::max(apex - left, min_left_right_count);
    int right_count = std::max(right - apex, min_left_right_count);
    max_left_count = std::max(max_left_count, left_count);
    max_right_count = std::max(max_right_count, right_count);
    entry.left_count = max_left_count;
    entry.right_count = max_right_count;
    double inte_sum = 0;
    double mass_sum = 0;
    for (int k = 0; k < num; k++) {
      inte[k] = inte[k] / std::sqrt(total_pwr);
      inte_sum += inte[k];
      mass_sum += inte[k] * mz[k];
    }
    entry.inte = std::move(inte);
    entry.avg_mass_delta = mass_sum / inte_sum - mz[0];
    entry.snr_mul_factor = inte_sum * inte_sum;
    entries_.push_back(std::move(entry));
  }
  if (entries_.empty()) {
    LOG_ERROR("The averagine table is empty!");
    exit(EXIT_FAILURE);
  }
}

using AveraginePtr = std::shared_ptr<Averagine>;

// ----------------------------------------------------------------------------
// LogMzPeak (FLASHHelperClasses::LogMzPeak)
// ----------------------------------------------------------------------------

// A spectrum peak on the logarithmic m/z axis, with the charge and the isotope
// index it was assigned in a peak group.
struct LogMzPeak {
  double mz = 0;
  double intensity = 0;
  double log_mz = 0;
  int abs_charge = 0;
  int isotope_index = 0;

  LogMzPeak() = default;

  LogMzPeak(double mz_value, double inte)
      : mz(mz_value),
        intensity(inte),
        log_mz(std::log(mz_value - kProtonMass)) {}

  double getUnchargedMass() const {
    if (abs_charge == 0) {
      return 0;
    }
    return (mz - kProtonMass) * abs_charge;
  }

  bool operator<(const LogMzPeak& a) const {
    if (log_mz == a.log_mz) {
      return intensity < a.intensity;
    }
    return log_mz < a.log_mz;
  }
};

using LogMzPeakVec = std::vector<LogMzPeak>;

// ----------------------------------------------------------------------------
// Isotope cosine (SpectralDeconvolution::getCosine,
// getIsotopeCosineAndIsoOffset)
// ----------------------------------------------------------------------------

// Cosine between the observed per-isotope intensities a[a_start, a_end) and
// the (L2-normalised) theoretical distribution b shifted by offset.
double getCosine(const std::vector<double>& a, int a_start, int a_end,
                 const std::vector<double>& b, int offset, int min_iso_len) {
  double n = 0;
  double a_norm = 0;
  a_start = std::max(0, a_start);
  a_end = std::min(static_cast<int>(a.size()), a_end);
  if (a_end - a_start < min_iso_len) {
    return 0;
  }
  for (int j = a_start; j < a_end; j++) {
    int i = j - offset;
    a_norm += a[j] * a[j];
    if (i >= static_cast<int>(b.size()) || i < 0 || b[i] <= 0) {
      continue;
    }
    n += a[j] * b[i];
  }
  if (a_norm <= 0) {
    return 0;
  }
  return n / std::sqrt(a_norm);
}

// The best cosine over the monoisotopic offsets -apex/4-1 .. apex/4+1
// (window_width 0 restricts the search to offset 0); the offset, in isotopes,
// that gave it is returned through offset. iso_int_shift is the index of
// isotope 0 in per_isotope_intensities (1: the array starts at isotope -1).
double getIsotopeCosineAndIsoOffset(
    double mono_mass, const std::vector<double>& per_isotope_intensities,
    int& offset, const Averagine& avg, int iso_int_shift, int window_width) {
  offset = 0;
  if (static_cast<int>(per_isotope_intensities.size()) <
      kMinIsoSize + iso_int_shift) {
    return 0;
  }
  const std::vector<double>& iso = avg.getInte(mono_mass);
  int right = avg.getApexIndex(mono_mass) / 4 + 1;
  int left = right;
  right += iso_int_shift;
  left -= iso_int_shift;
  double max_cos = -1000;
  int max_isotope_index = per_isotope_intensities.size();  // exclusive
  int min_isotope_index = -1;                              // inclusive
  for (int i = 0; i < max_isotope_index; i++) {
    if (per_isotope_intensities[i] <= 0) {
      continue;
    }
    if (min_isotope_index < 0) {
      min_isotope_index = i;
    }
  }
  if (max_isotope_index - min_isotope_index < kMinIsoSize) {
    return 0;
  }
  std::vector<std::pair<int, double>> offset_cos;
  offset_cos.reserve(right + left + 1);
  for (int tmp_offset = -left; tmp_offset <= right; tmp_offset++) {
    if (window_width >= 0 &&
        std::abs(tmp_offset - iso_int_shift) > window_width) {
      continue;
    }
    double tmp_cos = getCosine(per_isotope_intensities, min_isotope_index,
                               max_isotope_index, iso, tmp_offset, kMinIsoSize);
    offset_cos.emplace_back(tmp_offset, tmp_cos);
  }
  if (offset_cos.empty()) {
    return max_cos;
  }
  std::stable_sort(
      offset_cos.begin(), offset_cos.end(),
      [](const std::pair<int, double>& p1, const std::pair<int, double>& p2) {
        return p1.second > p2.second;
      });
  offset = offset_cos[0].first;
  max_cos = offset_cos[0].second;
  max_cos = std::max(max_cos, 0.0);
  offset -= iso_int_shift;
  return max_cos;
}

// The peak of a m/z-sorted list nearest to mz (MSSpectrum::findNearest).
std::size_t findNearest(const PeakPtrVec& peaks, double mz) {
  auto it = std::lower_bound(
      peaks.begin(), peaks.end(), mz,
      [](const PeakPtr& p, double value) { return p->getPosition() < value; });
  if (it == peaks.end()) {
    return peaks.size() - 1;
  }
  if (it == peaks.begin()) {
    return 0;
  }
  std::size_t idx = it - peaks.begin();
  if (mz - peaks[idx - 1]->getPosition() < peaks[idx]->getPosition() - mz) {
    return idx - 1;
  }
  return idx;
}

// ----------------------------------------------------------------------------
// PeakGroup (OpenMS PeakGroup and PeakGroupScoring)
// ----------------------------------------------------------------------------

// The peaks of one deconvoluted mass over its charge states and isotopes,
// with the monoisotopic mass, the intensity and the scores.
class PeakGroup {
 public:
  PeakGroup(int min_abs_charge, int max_abs_charge)
      : min_abs_charge_(min_abs_charge), max_abs_charge_(max_abs_charge) {}

  // Collect the raw peaks of the spectrum that fit mono_mass at the charges of
  // the group (PeakGroup::recruitAllPeaksInSpectrum with renew_signal_peaks):
  // the peaks within 0.8 x the tolerance of an isotope position become the
  // group's peaks, the others with a non-negative isotope index the noisy
  // peaks that are returned.
  LogMzPeakVec recruitAllPeaksInSpectrum(const PeakPtrVec& peaks, double tol,
                                         const Averagine& avg,
                                         double mono_mass);

  // Score the group for the current peaks (PeakGroup::updateQscore): the
  // charge range, the charge distribution score, the monoisotopic mass, the
  // isotope cosine with its offset, the per-charge cosines, the mass errors,
  // the SNRs and the QScore. Returns the isotope offset of the best cosine
  // (0 on the last pass, when only offset 0 is tried); the QScore is left 0
  // when the charge score, the cosine or the offset stop the scoring.
  int updateQscore(const LogMzPeakVec& noisy_peaks, const Averagine& avg,
                   double min_cos, double tol, bool is_low_charge,
                   bool is_last);

  // The intensity-weighted monoisotopic mass of the peaks and the per-isotope
  // intensities (PeakGroup::updateMonoMassAndIsotopeIntensities).
  void updateMonoMassAndIsotopeIntensities(double tol);

  void push_back(const LogMzPeak& p) { peaks_.push_back(p); }
  bool empty() const { return peaks_.empty(); }
  std::size_t size() const { return peaks_.size(); }
  const LogMzPeakVec& getPeaks() const { return peaks_; }
  void setPeaks(LogMzPeakVec peaks) { peaks_ = std::move(peaks); }
  void sort() { std::sort(peaks_.begin(), peaks_.end()); }

  void setMonoMass(double mass) { mono_mass_ = mass; }
  double getMonoMass() const { return mono_mass_; }
  double getIntensity() const { return intensity_; }
  int getMinAbsCharge() const { return min_abs_charge_; }
  int getMaxAbsCharge() const { return max_abs_charge_; }
  int getRepAbsCharge() const { return max_snr_abs_charge_; }
  double getIsotopeCosine() const { return isotope_cosine_score_; }
  double getChargeScore() const { return charge_score_; }
  double getSnr() const { return snr_; }
  double getQscore() const { return qscore_; }
  double getAvgPpmError() const { return avg_ppm_error_; }

  double getChargeSnr(int abs_charge) const {
    if (abs_charge < 0 ||
        abs_charge >= static_cast<int>(per_charge_snr_.size())) {
      return 0;
    }
    return per_charge_snr_[abs_charge];
  }

  double getChargeIsotopeCosine(int abs_charge) const {
    if (abs_charge < 0 ||
        abs_charge >= static_cast<int>(per_charge_cos_.size())) {
      return 0;
    }
    return per_charge_cos_[abs_charge];
  }

 private:
  LogMzPeakVec peaks_;
  LogMzPeakVec negative_iso_peaks_;
  std::vector<double> per_charge_sum_signal_squared_;
  std::vector<double> per_charge_noise_pwr_;
  std::vector<double> per_charge_cos_;
  std::vector<double> per_charge_int_;
  std::vector<double> per_charge_snr_;
  std::vector<double> per_isotope_int_;
  int min_abs_charge_ = 0;
  int max_abs_charge_ = -1;
  double mono_mass_ = -1.0;
  double intensity_ = 0;
  int max_snr_abs_charge_ = -1;
  double isotope_cosine_score_ = 0;
  double charge_score_ = 0;
  double qscore_ = 0;
  double avg_ppm_error_ = 0;
  double snr_ = 0;

  void clear() {
    peaks_.clear();
    negative_iso_peaks_.clear();
    intensity_ = 0;
  }

  void getPerIsotopeIntensities(std::vector<double>& intensities,
                                int& min_isotope_index, int& max_isotope_index,
                                int abs_charge, int min_negative_isotope_index,
                                double tol) const;
  void updatePerChargeInformation(const LogMzPeakVec& noisy_peaks, double tol,
                                  bool is_last);
  double getNoisePeakPower(const LogMzPeakVec& noisy_peaks, int z,
                           double tol) const;
  void updateChargeRange();
  void updateChargeFitScoreAndChargeIntensities(bool is_low_charge);
  void updatePerChargeCos(const Averagine& avg, double tol);
  double getPpmError(const LogMzPeak& p) const;
  void updateAvgMassError();
  void updateSnr(double mul_factor);
  double compQscore() const;
};

using PeakGroupVec = std::vector<PeakGroup>;

LogMzPeakVec PeakGroup::recruitAllPeaksInSpectrum(const PeakPtrVec& peaks,
                                                  double tol,
                                                  const Averagine& avg,
                                                  double mono_mass) {
  // not all peaks within the tolerance are signal peaks
  const double mul_tol = 0.8;
  LogMzPeakVec noisy_peaks;
  clear();
  if (mono_mass < 0) {
    return noisy_peaks;
  }
  // a charge range that is too narrow for the charge
  if (max_abs_charge_ - min_abs_charge_ < max_abs_charge_ / 20) {
    return noisy_peaks;
  }
  mono_mass_ = mono_mass;
  int max_isotope = avg.getLastIndex(mono_mass);
  int min_isotope = avg.getApexIndex(mono_mass) -
                    avg.getLeftCountFromApex(mono_mass) +
                    kMinNegativeIsotopeIndex;
  int max_signal_isotope = 0;
  min_isotope = std::max(kMinNegativeIsotopeIndex, min_isotope);
  int charge_num = max_abs_charge_ - min_abs_charge_ + 1;
  peaks_.reserve(max_isotope * charge_num * 2);
  noisy_peaks.reserve(max_isotope * charge_num * 2);
  int peak_num = peaks.size();
  // from the largest to the smallest charge
  for (int c = max_abs_charge_; c >= min_abs_charge_; c--) {
    if (c <= 0) {
      break;
    }
    double cmz = mono_mass / c + kProtonMass;
    double left_mz =
        (mono_mass - (1 - kMinNegativeIsotopeIndex) * kIsoDaDistance) / c +
        kProtonMass;
    int index = findNearest(peaks, left_mz * (1 - tol * mul_tol));
    double iso_delta = kIsoDaDistance / c;
    for (; index < peak_num; index++) {
      double pint = peaks[index]->getIntensity();
      if (pint <= 0) {
        continue;
      }
      double pmz = peaks[index]->getPosition();
      int iso_index = static_cast<int>(std::round((pmz - cmz) / iso_delta));
      if (iso_index > max_isotope) {
        break;
      }
      if (iso_index < min_isotope) {
        continue;
      }
      if (std::abs(pmz - cmz - iso_index * iso_delta) <= pmz * tol * mul_tol) {
        LogMzPeak p(pmz, pint);
        p.isotope_index = iso_index;
        p.abs_charge = c;
        if (iso_index < 0) {
          negative_iso_peaks_.push_back(p);
        } else {
          max_signal_isotope = std::max(max_signal_isotope, iso_index);
          peaks_.push_back(p);
        }
      } else if (iso_index >= 0) {
        LogMzPeak p(pmz, pint);
        p.isotope_index = static_cast<int>(std::floor((pmz - cmz) / iso_delta));
        p.abs_charge = c;
        noisy_peaks.push_back(p);
      }
    }
    if (index >= peak_num) {
      break;
    }
  }
  LogMzPeakVec kept_noisy_peaks;
  kept_noisy_peaks.reserve(noisy_peaks.size());
  for (const LogMzPeak& p : noisy_peaks) {
    if (p.isotope_index > max_signal_isotope) {
      continue;
    }
    kept_noisy_peaks.push_back(p);
  }
  return kept_noisy_peaks;
}

// PeakGroup::getPerIsotopeIntensities_: the summed intensity per isotope index
// (of one charge, or of all charges when abs_charge is 0), stored from
// isotope min_negative_isotope_index on, and smeared over the neighbouring
// isotopes when the tolerance in Da spans several of them (large masses only).
void PeakGroup::getPerIsotopeIntensities(std::vector<double>& intensities,
                                         int& min_isotope_index,
                                         int& max_isotope_index, int abs_charge,
                                         int min_negative_isotope_index,
                                         double tol) const {
  int da_tol = static_cast<int>(tol * peaks_[0].getUnchargedMass() / 2.0);
  min_isotope_index = 100000;
  max_isotope_index = -1;  // inclusive
  for (const LogMzPeak& p : peaks_) {
    max_isotope_index = std::max(max_isotope_index, p.isotope_index);
    min_isotope_index = std::min(min_isotope_index, p.isotope_index);
  }
  intensities.assign(max_isotope_index + 1 + da_tol - kMinNegativeIsotopeIndex,
                     0.0);
  double denom = 2.0 * std::pow(da_tol / 2.355, 2.0) / kIsoDaDistance;
  auto add_peak = [&](const LogMzPeak& peak) {
    if (abs_charge > 0 && peak.abs_charge != abs_charge) {
      return;
    }
    if (peak.isotope_index - min_negative_isotope_index >=
        static_cast<int>(intensities.size())) {
      return;
    }
    for (int margin = -da_tol; margin <= da_tol; margin++) {
      int index = peak.isotope_index + margin - min_negative_isotope_index;
      if (index < 0 || index >= static_cast<int>(intensities.size())) {
        continue;
      }
      intensities[index] +=
          denom > 0 ? peak.intensity * std::exp(-margin * margin / denom)
                    : peak.intensity;
    }
  };
  for (const LogMzPeak& peak : peaks_) {
    if (peak.isotope_index < 0) {
      continue;
    }
    add_peak(peak);
  }
  if (min_negative_isotope_index != 0) {
    for (const LogMzPeak& peak : negative_iso_peaks_) {
      if (peak.isotope_index - min_negative_isotope_index < 0) {
        continue;
      }
      add_peak(peak);
    }
  }
}

void PeakGroup::updateMonoMassAndIsotopeIntensities(double tol) {
  if (peaks_.empty()) {
    return;
  }
  int min_isotope_index;
  int max_isotope_index;
  sort();
  intensity_ = 0;
  double nominator = 0;
  getPerIsotopeIntensities(per_isotope_int_, min_isotope_index,
                           max_isotope_index, 0, kMinNegativeIsotopeIndex, tol);
  for (const LogMzPeak& p : peaks_) {
    if (p.isotope_index - kMinNegativeIsotopeIndex < 0) {
      continue;
    }
    nominator +=
        p.intensity * (p.getUnchargedMass() - p.isotope_index * kIsoDaDistance);
    intensity_ += p.intensity;
  }
  mono_mass_ = nominator / intensity_;
}

// PeakGroup::updatePerChargeInformation_: the intensity and the signal power
// (sum of squared per-isotope intensities) per charge, and the noise power
// per charge: on the last pass the path-based estimate of getNoisePeakPower,
// before it the plain sum of squared noisy peak intensities.
void PeakGroup::updatePerChargeInformation(const LogMzPeakVec& noisy_peaks,
                                           double tol, bool is_last) {
  per_charge_sum_signal_squared_.assign(1 + max_abs_charge_, 0.0);
  per_charge_int_.assign(1 + max_abs_charge_, 0.0);
  int max_iso = 0;
  for (const LogMzPeak& p : peaks_) {
    per_charge_int_[p.abs_charge] += p.intensity;
    max_iso = std::max(max_iso, p.isotope_index);
  }
  std::vector<std::vector<double>> per_charge_isotope_int(
      1 + max_abs_charge_, std::vector<double>(1 + max_iso, 0.0));
  for (const LogMzPeak& p : peaks_) {
    per_charge_isotope_int[p.abs_charge][p.isotope_index] += p.intensity;
  }
  for (int z = min_abs_charge_; z <= max_abs_charge_; z++) {
    for (int i = 0; i <= max_iso; i++) {
      double v = per_charge_isotope_int[z][i];
      per_charge_sum_signal_squared_[0] += v * v;
      per_charge_sum_signal_squared_[z] += v * v;
    }
  }
  per_charge_noise_pwr_.assign(1 + max_abs_charge_, 0.0);
  if (is_last) {
    for (int z = min_abs_charge_; z <= max_abs_charge_; z++) {
      per_charge_noise_pwr_[z] = getNoisePeakPower(noisy_peaks, z, tol);
    }
    per_charge_noise_pwr_[0] = getNoisePeakPower(noisy_peaks, 0, tol);
  } else {
    for (const LogMzPeak& p : noisy_peaks) {
      double pwr = p.intensity * p.intensity;
      if (p.abs_charge < min_abs_charge_ || p.abs_charge > max_abs_charge_) {
        continue;
      }
      per_charge_noise_pwr_[0] += pwr;
      per_charge_noise_pwr_[p.abs_charge] += pwr;
    }
  }
}

// PeakGroup::getNoisePeakPower_: the noise power of charge z (0: all
// charges). The noisy and the signal peaks are placed on the uncharged mass
// axis; series of peaks with one constant mass spacing (the spacing binned
// as bins 1..19 of normalised_distance / d * 12, d = 1, 2, 3, so that a
// series may skip one or two members) are the paths, and the strongest path
// per spacing is taken first, its unused peaks summed and the sum squared
// added to the power. The noisy peaks that no path took are grouped within
// the tolerance and add the square of each group's sum.
double PeakGroup::getNoisePeakPower(const LogMzPeakVec& noisy_peaks, int z,
                                    double tol) const {
  if (noisy_peaks.empty()) {
    return 0;
  }
  const std::size_t max_noisy_peak_number =
      z == 0 ? std::min(200, 40 * (max_abs_charge_ - min_abs_charge_ + 1)) : 40;
  const std::size_t bin_number_margin = 8;
  const std::size_t max_bin_number = bin_number_margin + 12;
  double threshold = -1;
  // (uncharged mass, intensity, is signal)
  struct MassPeak {
    double mass;
    double intensity;
    bool signal;
  };
  std::vector<MassPeak> all_peaks;
  all_peaks.reserve(max_noisy_peak_number + peaks_.size());
  std::size_t noise_peak_count = 0;
  for (const LogMzPeak& p : noisy_peaks) {
    if (z == 0 || p.abs_charge == z) {
      noise_peak_count++;
    }
  }
  if (noise_peak_count == 0) {
    return 0;
  }
  if (noise_peak_count > max_noisy_peak_number) {
    std::vector<double> intensities;
    intensities.reserve(noise_peak_count);
    for (const LogMzPeak& p : noisy_peaks) {
      if (z > 0 && p.abs_charge != z) {
        continue;
      }
      intensities.push_back(p.intensity);
    }
    std::sort(intensities.rbegin(), intensities.rend());
    threshold = intensities[max_noisy_peak_number];
  }
  for (const LogMzPeak& p : noisy_peaks) {
    if ((z > 0 && p.abs_charge != z) || p.intensity < threshold) {
      continue;
    }
    all_peaks.push_back({p.getUnchargedMass(), p.intensity, false});
  }
  if (all_peaks.empty()) {
    return 0;
  }
  for (const LogMzPeak& p : peaks_) {
    if ((z > 0 && p.abs_charge != z) || p.intensity < threshold) {
      continue;
    }
    all_peaks.push_back({p.getUnchargedMass(), p.intensity, true});
  }
  std::sort(
      all_peaks.begin(), all_peaks.end(),
      [](const MassPeak& a, const MassPeak& b) { return a.mass < b.mass; });
  double charge_noise_pwr = 0;
  std::size_t peak_num = all_peaks.size();
  std::vector<std::vector<std::size_t>> per_bin_edges(
      max_bin_number, std::vector<std::size_t>(peak_num, 0));
  // -2: the bin is empty, -1: used, >= 0: the start peak of the best path
  std::vector<int> per_bin_start_index(max_bin_number, -2);
  std::map<double, std::size_t> max_intensity_sum_to_bin;
  const std::vector<double> div_factors = {1.0, 2.0, 3.0};
  // the edges: a mass difference between two peaks, binned
  for (std::size_t i = 0; i < peak_num; i++) {
    double p1_mass = all_peaks[i].mass;
    std::vector<double> per_bin_error(max_bin_number, -1.0);
    for (std::size_t j = i + 1; j < peak_num; j++) {
      double normalized_dist = (all_peaks[j].mass - p1_mass) / kIsoDaDistance;
      // two signal peaks more than 0.75 isotopes apart are not connected
      if (all_peaks[i].signal && all_peaks[j].signal &&
          normalized_dist >= 0.75) {
        continue;
      }
      for (double d : div_factors) {
        double distance =
            normalized_dist / d * (max_bin_number - bin_number_margin);
        std::size_t bin = static_cast<std::size_t>(std::round(distance));
        if (bin == 0) {
          continue;
        }
        if (bin >= max_bin_number) {
          break;
        }
        per_bin_start_index[bin] = -1;
        double current_error =
            d * max_bin_number + std::abs(static_cast<double>(bin) - distance);
        if (per_bin_error[bin] >= 0 && per_bin_error[bin] < current_error) {
          continue;
        }
        per_bin_edges[bin][i] = j;
        per_bin_error[bin] = current_error;
      }
    }
  }
  // per bin, the path of the highest summed intensity
  for (std::size_t k = 0; k < max_bin_number; k++) {
    if (per_bin_start_index[k] == -2) {
      continue;
    }
    const std::vector<std::size_t>& edges = per_bin_edges[k];
    double max_sum_intensity = 0;
    for (std::size_t i = 0; i < edges.size(); i++) {
      if (edges[i] == 0) {
        break;
      }
      double sum_intensity = all_peaks[i].intensity;
      for (std::size_t j = edges[i]; j < edges.size(); j = edges[j]) {
        sum_intensity += all_peaks[j].intensity;
        if (j == 0) {
          break;
        }
      }
      if (max_sum_intensity < sum_intensity) {
        max_sum_intensity = sum_intensity;
        per_bin_start_index[k] = i;
      }
    }
    max_intensity_sum_to_bin[max_sum_intensity] = k;
  }
  std::vector<bool> unused(peak_num, true);
  // from the strongest path down, sum the intensities of the peaks no earlier
  // path used
  for (auto it = max_intensity_sum_to_bin.rbegin();
       it != max_intensity_sum_to_bin.rend(); ++it) {
    std::size_t bin = it->second;
    int index = per_bin_start_index[bin];
    if (index < 0) {
      continue;
    }
    const std::vector<std::size_t>& edges = per_bin_edges[bin];
    double ori_mass = all_peaks[index].mass;
    int ori_index = index;
    double sum_intensity = 0;
    auto follow_path = [&](int start) -> bool {
      std::size_t j = edges[start];
      if (j == 0) {
        return false;
      }
      if (unused[start]) {
        sum_intensity += all_peaks[start].intensity;
        unused[start] = false;
      } else {
        return false;
      }
      for (; j < edges.size(); j = edges[j]) {
        if (j == 0) {
          break;
        }
        if (unused[j]) {
          sum_intensity += all_peaks[j].intensity;
          unused[j] = false;
        } else {
          break;
        }
      }
      return true;
    };
    while (index < static_cast<int>(peak_num)) {
      double mass = all_peaks[index].mass;
      if (mass - ori_mass > tol / 2.0 * mass) {
        break;
      }
      if (!follow_path(index)) {
        break;
      }
      index++;
    }
    index = ori_index - 1;
    while (index >= 0) {
      double mass = all_peaks[index].mass;
      if (ori_mass - mass > tol / 2.0 * ori_mass) {
        break;
      }
      if (!follow_path(index)) {
        break;
      }
      index--;
    }
    charge_noise_pwr += sum_intensity * sum_intensity;
  }
  // the remaining noisy peaks, grouped within the tolerance
  double prev_mass = 0;
  double tmp_int_sum = 0;
  for (std::size_t index = 0; index < peak_num; index++) {
    if (!unused[index]) {
      continue;
    }
    const MassPeak& p = all_peaks[index];
    if (p.signal) {
      continue;
    }
    if (p.mass - prev_mass > p.mass * tol) {
      charge_noise_pwr += tmp_int_sum * tmp_int_sum;
      prev_mass = p.mass;
      tmp_int_sum = 0;
    }
    tmp_int_sum += p.intensity;
  }
  charge_noise_pwr += tmp_int_sum * tmp_int_sum;
  return charge_noise_pwr;
}

// PeakGroup::updateChargeRange_: from the charge of the highest signal-to-
// noise, extend the range up and down while the per-charge signal power over
// the noise power stays above a tenth of the maximum (at most 1), and drop
// the peaks outside the new range.
void PeakGroup::updateChargeRange() {
  int max_sig_charge = 0;
  double max_sig = 0;
  for (int z = min_abs_charge_; z <= max_abs_charge_; z++) {
    double tmp_snr = per_charge_int_[z] * per_charge_int_[z] /
                     (1 + per_charge_noise_pwr_[z]);
    if (max_sig < tmp_snr) {
      max_sig = tmp_snr;
      max_sig_charge = z;
    }
  }
  int new_max_abs_charge = max_sig_charge;
  int new_min_abs_charge = max_sig_charge;
  double threshold = std::min(max_sig / 10, 1.0);
  for (int z = max_sig_charge; z <= max_abs_charge_; z++) {
    double per_charge_signal_power = per_charge_int_[z] * per_charge_int_[z];
    if (per_charge_signal_power / (1 + per_charge_noise_pwr_[z]) < threshold) {
      break;
    }
    new_max_abs_charge = z;
  }
  for (int z = max_sig_charge; z >= min_abs_charge_; z--) {
    double per_charge_signal_power = per_charge_int_[z] * per_charge_int_[z];
    if (per_charge_signal_power / (1 + per_charge_noise_pwr_[z]) < threshold) {
      break;
    }
    new_min_abs_charge = z;
  }
  if (max_abs_charge_ != new_max_abs_charge ||
      min_abs_charge_ != new_min_abs_charge) {
    LogMzPeakVec new_peaks;
    new_peaks.reserve(peaks_.size());
    for (const LogMzPeak& p : peaks_) {
      if (p.abs_charge < new_min_abs_charge ||
          p.abs_charge > new_max_abs_charge) {
        continue;
      }
      new_peaks.push_back(p);
    }
    peaks_.swap(new_peaks);
    max_abs_charge_ = new_max_abs_charge;
    min_abs_charge_ = new_min_abs_charge;
  }
  if (min_abs_charge_ > max_abs_charge_) {
    clear();
  } else {
    sort();
  }
}

// PeakGroup::updateChargeFitScoreAndChargeIntensities_: 1 minus the fraction
// of the intensity that rises when walking away from the strongest charge
// (plus, for high charges, the part of a drop steeper than 30% per charge).
void PeakGroup::updateChargeFitScoreAndChargeIntensities(bool is_low_charge) {
  if (max_abs_charge_ == min_abs_charge_) {
    charge_score_ = 1;
    return;
  }
  double max_per_charge_intensity = 0;
  double summed_intensity = 0;
  int max_index = -1;
  int first_index = -1;
  int last_index = -1;
  for (int c = min_abs_charge_; c <= max_abs_charge_; c++) {
    summed_intensity += per_charge_int_[c];
    if (per_charge_int_[c] > 0) {
      if (first_index < 0) {
        first_index = c;
      }
      last_index = c;
    }
    if (max_per_charge_intensity > per_charge_int_[c]) {
      continue;
    }
    max_per_charge_intensity = per_charge_int_[c];
    max_index = c;
  }
  if (max_index < 0) {
    charge_score_ = 0;
    return;
  }
  first_index = first_index < 0 ? 0 : first_index;
  double p = 0;
  const double factor = 0.3;
  for (int c = max_index; c < last_index; c++) {
    double diff = per_charge_int_[c + 1] - per_charge_int_[c];
    if (diff > 0) {
      p += diff;
    } else if (!is_low_charge && diff < -per_charge_int_[c] * factor) {
      p -= diff + per_charge_int_[c] * factor;
    }
  }
  for (int c = max_index; c > first_index; c--) {
    double diff = per_charge_int_[c - 1] - per_charge_int_[c];
    if (diff > 0) {
      p += diff;
    } else if (!is_low_charge && diff < -per_charge_int_[c] * factor) {
      p -= diff + per_charge_int_[c] * factor;
    }
  }
  charge_score_ = std::max(0.0, 1.0 - p / summed_intensity);
}

// PeakGroup::updatePerChargeCos_: the isotope cosine of every charge state
// (offset 0).
void PeakGroup::updatePerChargeCos(const Averagine& avg, double tol) {
  const std::vector<double>& iso_dist = avg.getInte(mono_mass_);
  per_charge_cos_.assign(1 + max_abs_charge_, 0.0);
  if (min_abs_charge_ == max_abs_charge_) {
    per_charge_cos_[min_abs_charge_] = isotope_cosine_score_;
    return;
  }
  std::vector<double> current_per_isotope_intensities;
  for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_;
       abs_charge++) {
    int min_isotope_index;
    int max_isotope_index;
    getPerIsotopeIntensities(current_per_isotope_intensities, min_isotope_index,
                             max_isotope_index, abs_charge, 0, tol);
    per_charge_cos_[abs_charge] =
        getCosine(current_per_isotope_intensities, min_isotope_index,
                  max_isotope_index, iso_dist, 0, kMinIsoSize);
  }
}

double PeakGroup::getPpmError(const LogMzPeak& p) const {
  double mass = mono_mass_ + p.isotope_index * kIsoDaDistance;
  return (mass / p.abs_charge + kProtonMass - p.mz) / p.mz * 1e6;
}

// PeakGroup::updateAvgMassError_: the mean, over the isotope indices, of the
// absolute value of the mean signed ppm error of the peaks of that isotope.
void PeakGroup::updateAvgMassError() {
  std::vector<double> errors;
  int i = 0;
  double i_cntr = 0;
  double i_error = 0;
  for (const LogMzPeak& p : peaks_) {
    if (i != p.isotope_index) {
      if (i_cntr > 0) {
        errors.push_back(i_error / i_cntr);
      }
      i = p.isotope_index;
      i_cntr = 0;
      i_error = 0;
    }
    i_cntr++;
    i_error += getPpmError(p);
  }
  if (i_cntr > 0) {
    errors.push_back(i_error / i_cntr);
  }
  avg_ppm_error_ = 0;
  for (double e : errors) {
    avg_ppm_error_ += std::abs(e);
  }
  if (!errors.empty()) {
    avg_ppm_error_ /= errors.size();
  }
}

// PeakGroup::updateSNR_: per charge, the explained signal power
// (signal power x cosine^2) over the noise power plus the unexplained signal
// power; the group's SNR sums the charges' signal amplitudes.
void PeakGroup::updateSnr(double mul_factor) {
  per_charge_snr_.assign(1 + max_abs_charge_, 0.0);
  double total_nom = 0;
  double total_denom = 0;
  for (int c = min_abs_charge_; c <= max_abs_charge_; c++) {
    if (static_cast<int>(per_charge_cos_.size()) <= c) {
      continue;
    }
    double cos_squared = per_charge_cos_[c] * per_charge_cos_[c];
    double sig_pwr = per_charge_sum_signal_squared_[c] * cos_squared;
    double nom = 1e-6 * mul_factor + mul_factor * sig_pwr;
    double denom = 1e-6 * mul_factor + per_charge_noise_pwr_[c] +
                   (1 - cos_squared) * per_charge_sum_signal_squared_[c];
    per_charge_snr_[c] = denom <= 0 ? 0 : nom / denom;
    total_denom += (1 - cos_squared) * per_charge_sum_signal_squared_[c];
    total_nom += std::sqrt(sig_pwr);
  }
  snr_ = mul_factor * (1e-6 + total_nom * total_nom) /
         (1e-6 * mul_factor + total_denom + per_charge_noise_pwr_[0]);
  per_charge_sum_signal_squared_.clear();
  per_charge_noise_pwr_.clear();
}

// PeakGroupScoring::getQscore: a logistic regression on the isotope cosine,
// its difference to the representative charge's cosine, the log SNR of the
// representative charge, its difference to the log SNR of the group, and the
// average ppm error, with FLASHDeconv's fixed weights.
double PeakGroup::compQscore() const {
  if (peaks_.empty()) {
    return 0;
  }
  int rep = max_snr_abs_charge_;
  double rep_log_snr = std::log2(1 + getChargeSnr(rep));
  std::vector<double> fv = {isotope_cosine_score_,
                            isotope_cosine_score_ - getChargeIsotopeCosine(rep),
                            rep_log_snr, rep_log_snr - std::log2(1 + snr_),
                            avg_ppm_error_};
  double score = kQscoreWeights.back() + 0.5;
  for (std::size_t i = 0; i + 1 < kQscoreWeights.size(); i++) {
    score += fv[i] * kQscoreWeights[i];
  }
  return 1.0 / (1.0 + std::exp(score));
}

int PeakGroup::updateQscore(const LogMzPeakVec& noisy_peaks,
                            const Averagine& avg, double min_cos, double tol,
                            bool is_low_charge, bool is_last) {
  qscore_ = 0;
  if (empty()) {
    return 0;
  }
  updatePerChargeInformation(noisy_peaks, tol, is_last);
  updateChargeRange();
  updateChargeFitScoreAndChargeIntensities(is_low_charge);
  if (charge_score_ < kMinChargeScore) {
    return 0;
  }
  updateMonoMassAndIsotopeIntensities(tol);
  if (per_isotope_int_.empty() || max_abs_charge_ < min_abs_charge_) {
    return 0;
  }
  int h_offset;
  int window_width = is_last ? 0 : -1;
  isotope_cosine_score_ =
      getIsotopeCosineAndIsoOffset(mono_mass_, per_isotope_int_, h_offset, avg,
                                   -kMinNegativeIsotopeIndex, window_width);
  if (!is_last && h_offset != 0) {
    return h_offset;
  }
  if (isotope_cosine_score_ < min_cos) {
    return 0;
  }
  updatePerChargeCos(avg, tol);
  updateAvgMassError();
  updateSnr(avg.getSnrMulFactor(mono_mass_));
  max_snr_abs_charge_ = min_abs_charge_;
  for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_;
       abs_charge++) {
    if (getChargeSnr(abs_charge) > getChargeSnr(max_snr_abs_charge_)) {
      max_snr_abs_charge_ = abs_charge;
    }
  }
  qscore_ = compQscore();
  return h_offset;
}

// ----------------------------------------------------------------------------
// The deconvolution of one spectrum (SpectralDeconvolution)
// ----------------------------------------------------------------------------

class FlashDeconvSingleSp {
 public:
  // peak_list must be sorted by m/z; tol is the ppm tolerance divided by 1e6
  FlashDeconvSingleSp(const PeakPtrVec& peak_list, double tol, int max_charge,
                      double max_mass, const AveraginePtr& avg_ptr);

  // The peak groups of the spectrum, sorted by monoisotopic mass
  PeakGroupVec deconv();

 private:
  PeakPtrVec peak_list_;
  double tol_;
  double bin_mul_factor_;
  int max_charge_;
  double max_mass_;
  double min_mass_ = kMinMass;
  AveraginePtr avg_ptr_;

  LogMzPeakVec log_mz_peaks_;
  // -log(z) per charge index (the "universal pattern") and, per harmonic
  // charge and charge index, -log(z - floor(hc / 2) / hc)
  std::vector<double> universal_pattern_;
  std::vector<std::vector<double>> harmonic_pattern_matrix_;
  // the bin offsets of the patterns: mass bin = m/z bin + offset
  std::vector<int> binned_universal_pattern_;
  std::vector<std::vector<int>> binned_harmonic_patterns_;
  double mass_bin_min_value_ = 0;
  double mz_bin_min_value_ = 0;
  std::vector<bool> binned_log_mz_peaks_;
  std::vector<bool> binned_log_masses_;
  PeakGroupVec peak_groups_;

  static double getBinValue(std::size_t bin, double min_value, double factor) {
    return min_value + static_cast<double>(bin) / factor;
  }

  static std::size_t getBinNumber(double value, double min_value,
                                  double factor) {
    if (value < min_value) {
      return 0;
    }
    return static_cast<std::size_t>(std::round((value - min_value) * factor));
  }

  void filterLowPeaks();
  void updateLogMzPeaks();
  void setFilters();
  void binLogMzPeaks(std::size_t bin_number,
                     std::vector<double>& binned_log_mz_peak_intensities);
  void updateCandidateMassBins(std::vector<double>& mass_intensities,
                               const std::vector<double>& mz_intensities);
  void filterMassBins(const std::vector<double>& mass_intensities,
                      std::vector<int>& min_charge_idx,
                      std::vector<int>& max_charge_idx);
  void getCandidatePeakGroups(const std::vector<int>& min_charge_idx,
                              const std::vector<int>& max_charge_idx);
  void generatePeakGroupsFromSpectrum();
  void scoreAndFilterPeakGroups();
  void removeChargeErrorPeakGroups();
  void removeHarmonicPeakGroups();
  void removeOverlappingPeakGroups(double tol);
  static void removeDuplicateMasses(PeakGroupVec& groups);
};

FlashDeconvSingleSp::FlashDeconvSingleSp(const PeakPtrVec& peak_list,
                                         double tol, int max_charge,
                                         double max_mass,
                                         const AveraginePtr& avg_ptr)
    : peak_list_(peak_list),
      tol_(tol),
      bin_mul_factor_(1.0 / tol * kTolDivFactor),
      max_charge_(max_charge),
      max_mass_(max_mass),
      avg_ptr_(avg_ptr) {}

// FLASHDeconvAlgorithm::filterLowPeaks_: at most kMaxPeakNum peaks, none
// below a thousandth of the strongest one.
void FlashDeconvSingleSp::filterLowPeaks() {
  if (peak_list_.empty()) {
    return;
  }
  std::vector<double> intensities;
  intensities.reserve(peak_list_.size());
  double max_inte = 0;
  for (const PeakPtr& p : peak_list_) {
    intensities.push_back(p->getIntensity());
    max_inte = std::max(max_inte, p->getIntensity());
  }
  double threshold = max_inte * kMinPeakInteRatio;
  if (peak_list_.size() > kMaxPeakNum) {
    std::sort(intensities.rbegin(), intensities.rend());
    threshold = std::max(threshold, intensities[kMaxPeakNum]);
  }
  PeakPtrVec kept;
  kept.reserve(peak_list_.size());
  for (const PeakPtr& p : peak_list_) {
    if (p->getIntensity() >= threshold && p->getIntensity() > 0) {
      kept.push_back(p);
    }
  }
  peak_list_.swap(kept);
}

// SpectralDeconvolution::updateLogMzPeaks_
void FlashDeconvSingleSp::updateLogMzPeaks() {
  log_mz_peaks_.clear();
  log_mz_peaks_.reserve(peak_list_.size());
  for (const PeakPtr& p : peak_list_) {
    if (p->getIntensity() <= 0) {
      continue;
    }
    log_mz_peaks_.emplace_back(p->getPosition(), p->getIntensity());
  }
}

// SpectralDeconvolution::setFilters_
void FlashDeconvSingleSp::setFilters() {
  universal_pattern_.clear();
  int charge_range = max_charge_;
  for (int i = 0; i < charge_range; i++) {
    universal_pattern_.push_back(-std::log(i + 1));
  }
  harmonic_pattern_matrix_.assign(kHarmonicCharges.size(),
                                  std::vector<double>(charge_range, 0.0));
  for (std::size_t k = 0; k < kHarmonicCharges.size(); k++) {
    int hc = kHarmonicCharges[k];
    int n = hc / 2;
    for (int i = 0; i < charge_range; i++) {
      double a = i > 0 ? std::exp(-universal_pattern_[i - 1]) : 0;
      double b = std::exp(-universal_pattern_[i]);
      harmonic_pattern_matrix_[k][i] = -std::log(b - (b - a) * n / hc);
    }
  }
}

// SpectralDeconvolution::binLogMzPeaks_
void FlashDeconvSingleSp::binLogMzPeaks(
    std::size_t bin_number,
    std::vector<double>& binned_log_mz_peak_intensities) {
  binned_log_mz_peaks_.assign(bin_number, false);
  for (const LogMzPeak& p : log_mz_peaks_) {
    std::size_t bi = getBinNumber(p.log_mz, mz_bin_min_value_, bin_mul_factor_);
    if (bi >= bin_number) {
      break;
    }
    binned_log_mz_peaks_[bi] = true;
    binned_log_mz_peak_intensities[bi] += p.intensity;
  }
}

// SpectralDeconvolution::updateCandidateMassBins_: the mass bins supported by
// the spectrum. The m/z bins are visited from high to low m/z so that, for a
// mass, the charges come from small to large. A charge supports a mass when
// the previous charge supported it with an intensity within a factor (10 for
// charges up to 10, 5 + 50 / z above) or, for charges up to 10, when an
// isotope peak is present; a peak at a harmonic position (a fraction of a
// charge away) with an intensity between the two supports vetoes it. A mass
// bin is set when at least two charges (or z / 2 of them) support it, or one
// low charge with an isotope peak does.
void FlashDeconvSingleSp::updateCandidateMassBins(
    std::vector<double>& mass_intensities,
    const std::vector<double>& mz_intensities) {
  std::vector<std::size_t> mz_bin_indices;
  for (std::size_t i = 0; i < binned_log_mz_peaks_.size(); i++) {
    if (binned_log_mz_peaks_[i]) {
      mz_bin_indices.push_back(i);
    }
  }
  std::size_t h_charge_size = kHarmonicCharges.size();
  long bin_end = binned_log_masses_.size();
  std::vector<int> support_peak_count(binned_log_masses_.size(), 0);
  std::vector<int> prev_charges(binned_log_masses_.size(), max_charge_ + 2);
  std::vector<double> prev_intensities(binned_log_masses_.size(), 0.0);
  mass_intensities.assign(binned_log_masses_.size(), 0.0);
  std::vector<double> sub_max_h_intensity(h_charge_size, 0.0);
  const Averagine& avg = *avg_ptr_;
  const int tol_div = static_cast<int>(kTolDivFactor);
  for (auto iter = mz_bin_indices.rbegin(); iter != mz_bin_indices.rend();
       ++iter) {
    std::size_t mz_bin_index = *iter;
    const double intensity = mz_intensities[mz_bin_index];
    const double log_mz =
        getBinValue(mz_bin_index, mz_bin_min_value_, bin_mul_factor_);
    const double mz = std::exp(log_mz);
    const double iso_div_by_mz = kIsoDaDistance / mz;
    for (int j = 0; j < max_charge_; j++) {
      const long mass_bin_index =
          static_cast<long>(mz_bin_index) + binned_universal_pattern_[j];
      if (mass_bin_index < 0) {
        continue;
      }
      if (mass_bin_index >= bin_end) {
        break;
      }
      int& spc = support_peak_count[mass_bin_index];
      const int abs_charge = j + 1;
      double& prev_intensity = prev_intensities[mass_bin_index];
      int& prev_charge = prev_charges[mass_bin_index];
      const bool charge_not_continuous =
          prev_charge - j != -1 && prev_charge <= max_charge_;
      bool pass_first_check = false;
      const double highest_factor = 10.0;
      const double factor = abs_charge <= kLowCharge
                                ? highest_factor
                                : highest_factor / 2 + highest_factor / 2 *
                                                           kLowCharge /
                                                           abs_charge;
      const double hfactor = factor / 2.0;
      double intensity_ratio =
          prev_intensity <= 0 ? factor + 1 : intensity / prev_intensity;
      intensity_ratio =
          intensity_ratio < 1 ? 1.0 / intensity_ratio : intensity_ratio;
      double support_peak_intensity = 0;
      std::fill(sub_max_h_intensity.begin(), sub_max_h_intensity.end(), 0.0);
      if (charge_not_continuous || intensity_ratio > factor) {
        spc = 0;
      } else {
        pass_first_check = true;
        if (spc == 0 && abs_charge > kLowCharge) {
          support_peak_intensity = prev_intensity;
        }
      }
      // low charges: an isotope peak must be present
      if (!pass_first_check && abs_charge <= kLowCharge) {
        int d_end = avg.getApexIndex(mz * abs_charge) > 0 ? -1 : 1;
        for (int d = 1; d >= d_end; d -= 2) {
          bool iso_exist = false;
          int next_iso_bin = 0;
          const int nib = getBinNumber(log_mz + d * iso_div_by_mz / abs_charge,
                                       mz_bin_min_value_, bin_mul_factor_);
          const int nibr =
              abs_charge > 1
                  ? static_cast<int>(getBinNumber(
                        log_mz + d * iso_div_by_mz / (abs_charge - 1),
                        mz_bin_min_value_, bin_mul_factor_))
                  : 0;
          const int nibl =
              getBinNumber(log_mz + d * iso_div_by_mz / (abs_charge + 1),
                           mz_bin_min_value_, bin_mul_factor_);
          // the charges are not distinguishable
          if (std::abs(nib - nibr) < kTolDivFactor ||
              std::abs(nib - nibl) < kTolDivFactor) {
            break;
          }
          for (int t = -1; t < 2; t++) {
            int nibt = nib + t;
            if (std::abs(nibt - static_cast<int>(mz_bin_index)) >=
                    kTolDivFactor &&
                nibt > 0 &&
                nibt < static_cast<int>(binned_log_mz_peaks_.size()) &&
                binned_log_mz_peaks_[nibt]) {
              iso_exist = true;
              pass_first_check = true;
              if (next_iso_bin == 0 ||
                  mz_intensities[next_iso_bin] < mz_intensities[nibt]) {
                next_iso_bin = nibt;
              }
            }
          }
          // a harmonic peak between the peak and its isotope peak
          if (iso_exist) {
            const double h_threshold = intensity + mz_intensities[next_iso_bin];
            for (std::size_t k = 0; k < h_charge_size; k++) {
              const int hc = kHarmonicCharges[k];
              int harmonic_cntr = 0;
              if (hc * abs_charge > max_charge_) {
                break;
              }
              const int hdiff =
                  static_cast<int>(std::round(static_cast<double>(
                      next_iso_bin - static_cast<int>(mz_bin_index)))) /
                  hc * (hc / 2);
              const int next_harmonic_iso_bin =
                  static_cast<int>(mz_bin_index) + hdiff;
              if (std::abs(next_harmonic_iso_bin -
                           static_cast<int>(mz_bin_index)) >= kTolDivFactor &&
                  next_harmonic_iso_bin >= 0 &&
                  next_harmonic_iso_bin <
                      static_cast<int>(binned_log_mz_peaks_.size()) &&
                  binned_log_mz_peaks_[next_harmonic_iso_bin] &&
                  mz_intensities[next_harmonic_iso_bin] > h_threshold / 2 &&
                  mz_intensities[next_harmonic_iso_bin] < h_threshold * 2) {
                harmonic_cntr++;
                sub_max_h_intensity[k] += mz_intensities[next_harmonic_iso_bin];
              }
              if (harmonic_cntr > 0) {
                pass_first_check = false;
              }
            }
          }
          if (pass_first_check) {
            support_peak_intensity += mz_intensities[next_iso_bin];
          }
        }
        pass_first_check = pass_first_check &&
                           *std::max_element(sub_max_h_intensity.begin(),
                                             sub_max_h_intensity.end()) <= 0;
      }
      if (pass_first_check) {
        if (prev_charge - j == -1) {
          // consecutive charges: check the harmonic positions
          double max_intensity = intensity;
          double min_intensity = prev_intensity;
          if (prev_intensity <= 0) {
            max_intensity = intensity;
            min_intensity = intensity;
          } else if (min_intensity > max_intensity) {
            std::swap(min_intensity, max_intensity);
          }
          const double high_threshold = max_intensity * hfactor;
          const double low_threshold = min_intensity / hfactor;
          bool is_harmonic = false;
          int min_dis = tol_div + 1;
          for (std::size_t k = 0; k < h_charge_size; k++) {
            if (kHarmonicCharges[k] * abs_charge > max_charge_) {
              break;
            }
            double harmonic_intensity = 0;
            for (int t = -tol_div; t <= tol_div; t++) {
              long hmz_bin_index =
                  mass_bin_index - binned_harmonic_patterns_[k][j] + t;
              if (hmz_bin_index > 0 &&
                  hmz_bin_index != static_cast<long>(mz_bin_index) &&
                  hmz_bin_index <
                      static_cast<long>(binned_log_mz_peaks_.size()) &&
                  binned_log_mz_peaks_[hmz_bin_index]) {
                double h_intensity = mz_intensities[hmz_bin_index];
                if (h_intensity > low_threshold &&
                    h_intensity < high_threshold) {
                  is_harmonic = true;
                  if (std::abs(t) < min_dis) {
                    harmonic_intensity = h_intensity;
                    min_dis = std::abs(t);
                  }
                }
              }
            }
            sub_max_h_intensity[k] += harmonic_intensity;
          }
          if (!is_harmonic) {
            mass_intensities[mass_bin_index] +=
                intensity + support_peak_intensity;
            if (!binned_log_masses_[mass_bin_index]) {
              spc++;
              if (spc >= kMinSupportPeakCount || spc >= abs_charge / 2) {
                binned_log_masses_[mass_bin_index] = true;
              }
            }
          } else {
            mass_intensities[mass_bin_index] -= *std::max_element(
                sub_max_h_intensity.begin(), sub_max_h_intensity.end());
            if (spc > 0) {
              spc--;
            }
          }
        } else if (abs_charge <= kLowCharge) {
          // a low charge with an isotope peak
          mass_intensities[mass_bin_index] +=
              intensity + support_peak_intensity;
          if (!binned_log_masses_[mass_bin_index]) {
            spc++;
            binned_log_masses_[mass_bin_index] = true;
          }
        }
      }
      prev_intensity = intensity;
      prev_charge = j;
    }
  }
}

// SpectralDeconvolution::filterMassBins_: every m/z bin keeps the two
// candidate masses of the highest summed intensity (one when the best is more
// than twice the second); the charge index ranges of the kept masses are
// accumulated into min_charge_idx / max_charge_idx.
void FlashDeconvSingleSp::filterMassBins(
    const std::vector<double>& mass_intensities,
    std::vector<int>& min_charge_idx, std::vector<int>& max_charge_idx) {
  std::size_t bin_size = binned_log_masses_.size();
  min_charge_idx.assign(bin_size, std::numeric_limits<int>::max());
  max_charge_idx.assign(bin_size, std::numeric_limits<int>::min());
  std::vector<bool> to_skip(bin_size);
  for (std::size_t i = 0; i < bin_size; i++) {
    to_skip[i] = !binned_log_masses_[i];
    binned_log_masses_[i] = false;
  }
  const int select_top_n = 2;
  std::vector<long> max_indices(select_top_n, -1);
  std::vector<int> max_intensity_charge_idx(select_top_n, -1);
  for (std::size_t mz_bin_index = 0; mz_bin_index < binned_log_mz_peaks_.size();
       mz_bin_index++) {
    if (!binned_log_mz_peaks_[mz_bin_index]) {
      continue;
    }
    std::fill(max_indices.begin(), max_indices.end(), -1);
    std::fill(max_intensity_charge_idx.begin(), max_intensity_charge_idx.end(),
              -1);
    double max_intensity = 0;
    for (int j = 0; j < max_charge_; j++) {
      long mass_bin_index =
          static_cast<long>(mz_bin_index) + binned_universal_pattern_[j];
      if (mass_bin_index < 0) {
        continue;
      }
      if (mass_bin_index >= static_cast<long>(bin_size)) {
        break;
      }
      if (to_skip[mass_bin_index]) {
        continue;
      }
      double t = mass_intensities[mass_bin_index];
      if (t == 0) {
        continue;
      }
      if (max_intensity == 0 || max_intensity < t) {
        bool big_diff = t > max_intensity * 2.0;
        for (int i = select_top_n - 1; i > 0; i--) {
          max_indices[i] = big_diff ? -1 : max_indices[i - 1];
          max_intensity_charge_idx[i] = max_intensity_charge_idx[i - 1];
        }
        max_indices[0] = mass_bin_index;
        max_intensity_charge_idx[0] = j;
        max_intensity = t;
      }
    }
    for (int i = 0; i < select_top_n; i++) {
      long max_index = max_indices[i];
      if (max_index >= 0 && max_index < static_cast<long>(bin_size)) {
        min_charge_idx[max_index] =
            std::min(min_charge_idx[max_index], max_intensity_charge_idx[i]);
        max_charge_idx[max_index] =
            std::max(max_charge_idx[max_index], max_intensity_charge_idx[i]);
        binned_log_masses_[max_index] = true;
      }
    }
  }
}

// SpectralDeconvolution::getCandidatePeakGroups_: for every candidate mass
// bin and every charge of its range, the strongest peak in the m/z bin (a
// local maximum among the neighbouring bins) seeds the collection of the
// isotope peaks to the right and left within the averagine span; the peaks
// are then re-indexed against the strongest peak of the group, which is set
// to the averagine apex, and a group with at least two distinct isotope
// indices becomes a candidate peak group.
void FlashDeconvSingleSp::getCandidatePeakGroups(
    const std::vector<int>& min_charge_idx,
    const std::vector<int>& max_charge_idx) {
  const Averagine& avg = *avg_ptr_;
  int charge_range = max_charge_;
  std::size_t mass_bin_size = binned_log_masses_.size();
  int log_mz_peak_size = log_mz_peaks_.size();
  std::vector<int> current_peak_index(charge_range, 0);
  std::vector<std::size_t> peak_bin_numbers(log_mz_peak_size);
  for (int i = 0; i < log_mz_peak_size; i++) {
    peak_bin_numbers[i] = getBinNumber(log_mz_peaks_[i].log_mz,
                                       mz_bin_min_value_, bin_mul_factor_);
  }
  for (std::size_t mass_bin_index = 0; mass_bin_index < mass_bin_size;
       mass_bin_index++) {
    if (!binned_log_masses_[mass_bin_index]) {
      continue;
    }
    double log_m =
        getBinValue(mass_bin_index, mass_bin_min_value_, bin_mul_factor_);
    double mass = std::exp(log_m);
    PeakGroup pg(1, max_charge_idx[mass_bin_index] + 1);
    std::size_t right_index = avg.getRightCountFromApex(mass);
    std::size_t left_index = avg.getLeftCountFromApex(mass);
    for (int j = min_charge_idx[mass_bin_index];
         j <= max_charge_idx[mass_bin_index]; j++) {
      int max_peak_index = -1;
      int abs_charge = j + 1;
      int bin_offset = binned_universal_pattern_[j];
      if (static_cast<long>(mass_bin_index) < bin_offset) {
        continue;
      }
      std::size_t b_index = mass_bin_index - bin_offset;  // m/z bin
      int& cpi = current_peak_index[j];
      double max_intensity = -1;
      while (cpi < log_mz_peak_size - 1) {
        if (peak_bin_numbers[cpi] == b_index) {
          double intensity = log_mz_peaks_[cpi].intensity;
          if (intensity > max_intensity) {
            max_intensity = intensity;
            max_peak_index = cpi;
          }
        } else if (peak_bin_numbers[cpi] > b_index) {
          break;
        }
        cpi++;
      }
      if (max_peak_index < 0) {
        continue;
      }
      // a local maximum
      if (max_peak_index > 0 &&
          peak_bin_numbers[max_peak_index - 1] == b_index - 1 &&
          log_mz_peaks_[max_peak_index - 1].intensity > max_intensity) {
        continue;
      }
      if (max_peak_index < log_mz_peak_size - 1 &&
          peak_bin_numbers[max_peak_index + 1] == b_index + 1 &&
          log_mz_peaks_[max_peak_index + 1].intensity > max_intensity) {
        continue;
      }
      const double mz = log_mz_peaks_[max_peak_index].mz;
      const double iso_delta = kIsoDaDistance / abs_charge;
      double mz_delta =
          std::min(kMaxMassDaltonTolerance / abs_charge, tol_ * mz);
      double max_mz = mz;
      for (int peak_index = max_peak_index; peak_index < log_mz_peak_size;
           peak_index++) {
        const double observed_mz = log_mz_peaks_[peak_index].mz;
        double mz_diff = observed_mz - mz;
        int tmp_i = static_cast<int>(std::round(mz_diff / iso_delta));
        if (observed_mz - max_mz > right_index * iso_delta + mz_delta) {
          break;
        }
        if (std::abs(mz_diff - tmp_i * iso_delta) < mz_delta) {
          std::size_t bin = peak_bin_numbers[peak_index] + bin_offset;
          if (bin < mass_bin_size) {
            LogMzPeak p(log_mz_peaks_[peak_index]);
            p.abs_charge = abs_charge;
            p.isotope_index = tmp_i;
            pg.push_back(p);
          }
        }
      }
      for (int peak_index = max_peak_index - 1; peak_index >= 0; peak_index--) {
        const double observed_mz = log_mz_peaks_[peak_index].mz;
        double mz_diff = observed_mz - mz;
        int tmp_i = static_cast<int>(std::round(mz_diff / iso_delta));
        if (max_mz - observed_mz > left_index * iso_delta + mz_delta) {
          break;
        }
        if (std::abs(mz_diff - tmp_i * iso_delta) < mz_delta) {
          std::size_t bin = peak_bin_numbers[peak_index] + bin_offset;
          if (bin < mass_bin_size) {
            LogMzPeak p(log_mz_peaks_[peak_index]);
            p.abs_charge = abs_charge;
            p.isotope_index = tmp_i;
            pg.push_back(p);
          }
        }
      }
    }
    if (pg.empty()) {
      continue;
    }
    double max_intensity = -1.0;
    double t_mass = 0;
    for (const LogMzPeak& p : pg.getPeaks()) {
      if (max_intensity < p.intensity) {
        max_intensity = p.intensity;
        t_mass = p.getUnchargedMass();
      }
    }
    double iso_tolerance = tol_ * t_mass;
    int min_off = 10000;
    int max_off = -1;
    int apex_index = avg.getApexIndex(t_mass);
    LogMzPeakVec new_peaks;
    new_peaks.reserve(pg.size());
    for (LogMzPeak p : pg.getPeaks()) {
      p.isotope_index = static_cast<int>(
          std::round((p.getUnchargedMass() - t_mass) / kIsoDaDistance));
      if (std::abs(t_mass - p.getUnchargedMass() +
                   kIsoDaDistance * p.isotope_index) > iso_tolerance) {
        continue;
      }
      p.isotope_index += apex_index;
      new_peaks.push_back(p);
      min_off = std::min(min_off, p.isotope_index);
      max_off = std::max(max_off, p.isotope_index);
    }
    if (min_off == max_off) {
      continue;
    }
    pg.setPeaks(std::move(new_peaks));
    pg.updateMonoMassAndIsotopeIntensities(tol_);
    if (pg.getMonoMass() < min_mass_ || pg.getMonoMass() > max_mass_) {
      continue;
    }
    pg.sort();
    peak_groups_.push_back(std::move(pg));
  }
}

// SpectralDeconvolution::generatePeakGroupsFromSpectrum_
void FlashDeconvSingleSp::generatePeakGroupsFromSpectrum() {
  peak_groups_.clear();
  const Averagine& avg = *avg_ptr_;
  int charge_range = max_charge_;
  int tmp_peak_cntr = std::max(0, max_charge_ - kMinSupportPeakCount);
  double mass_bin_max_value = std::min(
      log_mz_peaks_.back().log_mz - universal_pattern_[tmp_peak_cntr],
      std::log(max_mass_ + avg.getRightCountFromApex(max_mass_) + 1.0));
  mass_bin_min_value_ =
      std::log(std::max(1.0, 50 - avg.getAverageMassDelta(50)));
  mz_bin_min_value_ = log_mz_peaks_[0].log_mz;
  double mz_bin_max_value = log_mz_peaks_.back().log_mz;
  std::size_t mass_bin_number =
      getBinNumber(mass_bin_max_value, mass_bin_min_value_, bin_mul_factor_) +
      1;
  binned_universal_pattern_.clear();
  for (int i = 0; i < charge_range; i++) {
    binned_universal_pattern_.push_back(static_cast<int>(std::round(
        (mz_bin_min_value_ - universal_pattern_[i] - mass_bin_min_value_) *
        bin_mul_factor_)));
  }
  binned_harmonic_patterns_.assign(kHarmonicCharges.size(),
                                   std::vector<int>(charge_range, 0));
  for (std::size_t k = 0; k < kHarmonicCharges.size(); k++) {
    for (int i = 0; i < charge_range; i++) {
      binned_harmonic_patterns_[k][i] = static_cast<int>(
          std::round((mz_bin_min_value_ - harmonic_pattern_matrix_[k][i] -
                      mass_bin_min_value_) *
                     bin_mul_factor_));
    }
  }
  std::size_t mz_bin_number =
      getBinNumber(mz_bin_max_value, mz_bin_min_value_, bin_mul_factor_) + 1;
  std::vector<double> mz_intensities(mz_bin_number, 0.0);
  binLogMzPeaks(mz_bin_number, mz_intensities);
  binned_log_masses_.assign(mass_bin_number, false);
  std::vector<double> mass_intensities;
  updateCandidateMassBins(mass_intensities, mz_intensities);
  std::vector<int> min_charge_idx;
  std::vector<int> max_charge_idx;
  filterMassBins(mass_intensities, min_charge_idx, max_charge_idx);
  getCandidatePeakGroups(min_charge_idx, max_charge_idx);
  scoreAndFilterPeakGroups();
}

// keep one peak group per monoisotopic mass (the groups are sorted by mass)
void FlashDeconvSingleSp::removeDuplicateMasses(PeakGroupVec& groups) {
  PeakGroupVec filtered;
  filtered.reserve(groups.size());
  for (PeakGroup& pg : groups) {
    if (filtered.empty() || filtered.back().getMonoMass() < pg.getMonoMass()) {
      filtered.push_back(std::move(pg));
    }
  }
  groups.swap(filtered);
}

// SpectralDeconvolution::removeChargeErrorPeakGroups_: a peak group whose
// peaks are shared, at another charge, with peak groups of a higher
// per-charge SNR for at least half of its intensity is a charge error.
void FlashDeconvSingleSp::removeChargeErrorPeakGroups() {
  std::map<double, std::set<int>> peak_to_pgs;
  std::map<double, double> mz_to_intensities;
  std::vector<double> overlap_intensity(peak_groups_.size(), 0.0);
  for (std::size_t i = 0; i < peak_groups_.size(); i++) {
    for (const LogMzPeak& p : peak_groups_[i].getPeaks()) {
      peak_to_pgs[p.mz].insert(i);
      mz_to_intensities[p.mz] = p.intensity;
    }
  }
  const double high_snr = 10;
  for (const auto& e : peak_to_pgs) {
    const std::set<int>& pg_is = e.second;
    double pmz = e.first;
    double pint = mz_to_intensities[pmz];
    if (pg_is.size() == 1) {
      continue;
    }
    for (int i : pg_is) {
      bool is_overlap = false;
      double mass1 = peak_groups_[i].getMonoMass();
      int repz1 = static_cast<int>(std::round(mass1 / (pmz - kProtonMass)));
      // an almost perfect signal is not filtered
      if (peak_groups_[i].getSnr() > high_snr &&
          peak_groups_[i].getChargeSnr(repz1) > high_snr) {
        continue;
      }
      for (int j : pg_is) {
        if (i == j) {
          continue;
        }
        double mass2 = peak_groups_[j].getMonoMass();
        int repz2 = static_cast<int>(std::round(mass2 / (pmz - kProtonMass)));
        if (repz1 == repz2) {
          continue;
        }
        if (peak_groups_[i].getChargeSnr(repz1) >
            peak_groups_[j].getChargeSnr(repz2)) {
          continue;
        }
        is_overlap = true;
        break;
      }
      if (is_overlap) {
        overlap_intensity[i] += pint;
      }
    }
  }
  PeakGroupVec filtered;
  filtered.reserve(peak_groups_.size());
  const double mul_factor = 0.5;
  for (std::size_t i = 0; i < peak_groups_.size(); i++) {
    if (peak_groups_[i].getRepAbsCharge() > max_charge_) {
      continue;
    }
    if (overlap_intensity[i] >= peak_groups_[i].getIntensity() * mul_factor) {
      continue;
    }
    filtered.push_back(std::move(peak_groups_[i]));
  }
  peak_groups_.swap(filtered);
}

// The harmonic removal of scoreAndFilterPeakGroups_: a peak group whose mass
// times 2, 3, 5, 7 or 11 also deconvolutes, with a higher SNR, is the
// harmonic of that mass.
void FlashDeconvSingleSp::removeHarmonicPeakGroups() {
  const Averagine& avg = *avg_ptr_;
  PeakGroupVec filtered;
  filtered.reserve(peak_groups_.size());
  for (PeakGroup& peak_group : peak_groups_) {
    bool pass = true;
    int z1 = peak_group.getMinAbsCharge();
    int z2 = peak_group.getMaxAbsCharge();
    for (int hz : kHarmonicCharges) {
      PeakGroup pg(std::min(max_charge_, z1 * hz),
                   std::min(max_charge_, z2 * hz));
      // the harmonic charges are clipped too much
      if (z2 - z1 > pg.getMaxAbsCharge() - pg.getMinAbsCharge()) {
        break;
      }
      pg.setMonoMass(peak_group.getMonoMass() * hz);
      LogMzPeakVec nps =
          pg.recruitAllPeaksInSpectrum(peak_list_, tol_, avg, pg.getMonoMass());
      pg.updateQscore(nps, avg, kMinIsotopeCosine, tol_,
                      (z1 + z2) * hz < 2 * kLowCharge, true);
      if (pg.getQscore() > 0 && pg.getSnr() > peak_group.getSnr()) {
        pass = false;
        break;
      }
    }
    if (pass) {
      filtered.push_back(std::move(peak_group));
    }
  }
  peak_groups_.swap(filtered);
}

// SpectralDeconvolution::removeOverlappingPeakGroups_: of the peak groups
// whose masses agree within tol (a greedy clustering along the sorted
// masses), only the one of the highest SNR is kept.
void FlashDeconvSingleSp::removeOverlappingPeakGroups(double tol) {
  if (peak_groups_.empty()) {
    return;
  }
  std::vector<std::size_t> kept;
  double start_mass = peak_groups_[0].getMonoMass();
  double local_max_snr = -1.0;
  std::size_t local_max_index = 0;
  std::size_t last_local_max_index = peak_groups_.size();
  for (std::size_t i = 0; i < peak_groups_.size(); i++) {
    double mass = peak_groups_[i].getMonoMass();
    if (mass - start_mass > mass * tol) {
      if (last_local_max_index != local_max_index) {
        kept.push_back(local_max_index);
      }
      last_local_max_index = local_max_index;
      start_mass = mass;
      local_max_snr = -1.0;
    }
    double snr = peak_groups_[i].getSnr();
    if (local_max_snr < snr) {
      local_max_snr = snr;
      local_max_index = i;
    }
  }
  if (local_max_snr >= 0 && last_local_max_index != local_max_index) {
    kept.push_back(local_max_index);
  }
  PeakGroupVec filtered;
  filtered.reserve(kept.size());
  for (std::size_t i : kept) {
    filtered.push_back(std::move(peak_groups_[i]));
  }
  peak_groups_.swap(filtered);
}

// SpectralDeconvolution::scoreAndFilterPeakGroups_
void FlashDeconvSingleSp::scoreAndFilterPeakGroups() {
  const Averagine& avg = *avg_ptr_;
  PeakGroupVec scored;
  scored.reserve(peak_groups_.size());
  for (PeakGroup& peak_group : peak_groups_) {
    int offset = 0;
    bool mass_determined = false;
    double prev_cos = 0;
    // the isotope pattern is matched again until the cosine stops improving
    // and the monoisotopic mass settles
    for (int k = 0; k < kRefineIterationNum; k++) {
      LogMzPeakVec noisy_peaks = peak_group.recruitAllPeaksInSpectrum(
          peak_list_, tol_, avg,
          peak_group.getMonoMass() + offset * kIsoDaDistance);
      int z1 = peak_group.getMinAbsCharge();
      int z2 = peak_group.getMaxAbsCharge();
      bool is_low_charge = (z1 + z2) < 2 * kLowCharge;
      offset = peak_group.updateQscore(noisy_peaks, avg, kMinIsotopeCosine,
                                       tol_, is_low_charge, false);
      if (prev_cos > peak_group.getIsotopeCosine() || offset == 0 ||
          k >= kRefineIterationNum - 1) {
        if (peak_group.getChargeSnr(peak_group.getRepAbsCharge()) < kMinSnr) {
          break;
        }
        if (offset != 0) {
          noisy_peaks = peak_group.recruitAllPeaksInSpectrum(
              peak_list_, tol_, avg,
              peak_group.getMonoMass() + offset * kIsoDaDistance);
        }
        peak_group.updateQscore(noisy_peaks, avg, kMinIsotopeCosine, tol_,
                                is_low_charge, true);
        mass_determined = true;
        break;
      }
      prev_cos = peak_group.getIsotopeCosine();
    }
    if (!mass_determined || peak_group.empty() || peak_group.getQscore() <= 0 ||
        peak_group.getMonoMass() < min_mass_ ||
        peak_group.getMonoMass() > max_mass_) {
      continue;
    }
    int z1 = peak_group.getMinAbsCharge();
    int z2 = peak_group.getMaxAbsCharge();
    if (z1 > kLowCharge && (z2 - z1) < kMinSupportPeakCount) {
      continue;
    }
    scored.push_back(std::move(peak_group));
  }
  peak_groups_.swap(scored);
  if (peak_groups_.empty()) {
    return;
  }
  auto cmp_mass = [](const PeakGroup& a, const PeakGroup& b) {
    if (a.getMonoMass() == b.getMonoMass()) {
      return a.getIntensity() > b.getIntensity();
    }
    return a.getMonoMass() < b.getMonoMass();
  };
  // the duplicates are removed in the candidate order, then the groups are
  // sorted by mass, as in SpectralDeconvolution
  removeDuplicateMasses(peak_groups_);
  std::sort(peak_groups_.begin(), peak_groups_.end(), cmp_mass);
  removeChargeErrorPeakGroups();
  removeHarmonicPeakGroups();
  removeDuplicateMasses(peak_groups_);
  std::sort(peak_groups_.begin(), peak_groups_.end(), cmp_mass);
  removeOverlappingPeakGroups(tol_ * 1.2);
  PeakGroupVec filtered;
  filtered.reserve(peak_groups_.size());
  for (PeakGroup& pg : peak_groups_) {
    if (pg.getQscore() >= kMinQscore &&
        pg.getChargeSnr(pg.getRepAbsCharge()) > kMinSnr) {
      filtered.push_back(std::move(pg));
    }
  }
  peak_groups_.swap(filtered);
  std::sort(peak_groups_.begin(), peak_groups_.end(), cmp_mass);
}

PeakGroupVec FlashDeconvSingleSp::deconv() {
  peak_groups_.clear();
  filterLowPeaks();
  updateLogMzPeaks();
  if (log_mz_peaks_.empty() || max_charge_ < 1) {
    return peak_groups_;
  }
  setFilters();
  generatePeakGroupsFromSpectrum();
  return peak_groups_;
}

// ----------------------------------------------------------------------------
// From peak groups to TopFD envelopes
// ----------------------------------------------------------------------------

// A peak group becomes a MatchEnv of its representative (highest per-charge
// SNR) charge: the theoretical envelope of its monoisotopic mass from EnvBase,
// scaled so that its intensity sum is the peak group's intensity over all
// charges, matched against the spectrum peaks (ExpEnv). The QScore is stored
// as the EnvCNN score (the score the msalign file reports and sorts by) and
// the isotope cosine as the MsDeconv score.
//
// bound_inte is the spectrum's noise floor (0 when TopFD's own noise
// estimation is off); after scaling, the theoretical envelope is trimmed to
// the peaks whose intensity stays above it, exactly as the standard
// DeconvMs2Process pipeline trims its own final envelope
// (env_detect::compFullEnv, via getCandidateEnv: getSubEnv(1.0, bound_inte,
// INT_MAX, INT_MAX)). Skipping this step left every envelope carrying the
// averagine table's full isotope tail -- peaks orders of magnitude below the
// apex that no real spectrum could show above noise -- which inflated the
// envelopes this class writes to the SQLite database (toppic_test experiment,
// median 14 peaks per MS/MS envelope against DeconvMs2Process's 4 on the same
// data). min_inte is the matching threshold ExpEnv itself applies to observed
// peaks (DeconvData::initMinInte's min_inte_, equal to bound_inte whenever
// -ms2-sn-ratio is left at its default of 1).
MatchEnvPtrVec toMatchEnvs(const PeakGroupVec& peak_groups,
                           const PeakPtrVec& peak_list,
                           const EnvParaPtr& env_para_ptr, double min_inte,
                           double bound_inte) {
  MatchEnvPtrVec envs;
  envs.reserve(peak_groups.size());
  for (const PeakGroup& pg : peak_groups) {
    double mono_mass = pg.getMonoMass();
    int charge = pg.getRepAbsCharge();
    if (charge < 1) {
      continue;
    }
    EnvPtr theo_env = EnvBase::getEnvByMonoMass(mono_mass, charge);
    if (theo_env == nullptr) {
      continue;
    }
    double theo_sum = theo_env->compInteSum();
    if (theo_sum > 0) {
      theo_env->changeIntensity(pg.getIntensity() / theo_sum);
    }
    theo_env =
        theo_env->getSubEnv(1.0, bound_inte, std::numeric_limits<int>::max(),
                            std::numeric_limits<int>::max());
    double tolerance = env_para_ptr->getMzTolerance(charge);
    ExpEnvPtr exp_env =
        std::make_shared<ExpEnv>(peak_list, theo_env, tolerance, min_inte);
    int mass_group = env_para_ptr->getMassGroup(mono_mass);
    MatchEnvPtr match_env =
        std::make_shared<MatchEnv>(mass_group, theo_env, exp_env);
    match_env->setEnvcnnScore(pg.getQscore());
    match_env->setMsdeconvScore(pg.getIsotopeCosine());
    envs.push_back(match_env);
  }
  return envs;
}

// ----------------------------------------------------------------------------
// The processing of a file (the same as deconv_ms2_process)
// ----------------------------------------------------------------------------

// One deconvoluted MS/MS spectrum's data, buffered in memory so the slow SQLite
// writes can be deferred and done sequentially after the parallel
// deconvolution finishes.
struct Ms2SqlRecord {
  MzmlMsPtr ms_ptr;
  MatchEnvPtrVec deconv_envs;
};

// One record buffer per worker thread (indexed by writer id) so threads append
// without locking, plus an approximate byte counter shared across threads so
// the main thread can cap the buffer's memory use and flush it when full.
struct Ms2SqlBuffer {
  explicit Ms2SqlBuffer(int thread_num) : thread_records(thread_num) {}
  std::vector<std::vector<Ms2SqlRecord>> thread_records;
  std::atomic<std::size_t> bytes{0};
};
using Ms2SqlBufferPtr = std::shared_ptr<Ms2SqlBuffer>;

// Approximate bytes held per buffered peak, used only to bound the buffer.
constexpr std::size_t BYTES_PER_PEAK = 64;

// Hard ceiling for the in-memory SQL buffer.
constexpr std::size_t MAX_SQL_BUFFER_BYTES =
    4ULL * 1024 * 1024 * 1024;  // 4 GiB

// Start draining at 90% of the ceiling, leaving headroom for the records that
// in-flight worker tasks keep appending while the buffer is being drained.
constexpr std::size_t SQL_BUFFER_FLUSH_BYTES = MAX_SQL_BUFFER_BYTES * 9 / 10;

std::size_t estimateRecordBytes(const Ms2SqlRecord& record) {
  std::size_t peak_num = record.ms_ptr->size();
  for (const auto& env : record.deconv_envs) {
    peak_num += env->getExpEnvPtr()->getPeakNum();
  }
  return peak_num * BYTES_PER_PEAK;
}

std::string updateMsTwoMsg(const MsHeaderPtr& header_ptr, int scan_cnt,
                           int total_scan_num) {
  std::string percentage = std::to_string(scan_cnt * 100 / total_scan_num);
  std::string msg = "Processing MS/MS spectrum scan " +
                    std::to_string(header_ptr->getFirstScanNum()) + " ...";
  while (msg.length() < 40) {
    msg += " ";
  }
  msg = msg + percentage + "% finished.";
  return msg;
}

void deconvMsTwo(const MzmlMsPtr& ms_ptr, SpecFeaturePtrVec sp_feat_ptr_vec,
                 const TopfdParaPtr& topfd_para_ptr,
                 const AveraginePtr& avg_ptr, const EnvParaPtr& env_para_ptr,
                 const MsAlignWriterPtrVec& ms2_writer_ptr_vec,
                 const SimpleThreadPoolPtr& pool_ptr,
                 const Ms2SqlBufferPtr& sql_buffer_ptr) {
  // 1. Find max_mass and max_charge: the precursor's, capped by the arguments
  // (FLASHDeconv caps the MSn mass and charge by its precursor peak group)
  double max_mass = 0;
  int max_charge = 1;
  for (std::size_t i = 0; i < sp_feat_ptr_vec.size(); i++) {
    double mass = sp_feat_ptr_vec[i]->getPrecMonoMass();
    if (mass > max_mass) {
      max_mass = mass;
    }
    int charge = sp_feat_ptr_vec[i]->getPrecCharge();
    if (charge > max_charge) {
      max_charge = charge;
    }
  }
  double arg_max_mass = topfd_para_ptr->getMaxMass();
  if (arg_max_mass < max_mass) {
    max_mass = arg_max_mass;
  }
  int arg_max_charge = topfd_para_ptr->getMaxCharge();
  if (arg_max_charge < max_charge) {
    max_charge = arg_max_charge;
  }
  // 2. Deconvolute the spectrum
  PeakPtrVec peak_list = ms_ptr->getPeakPtrVec();
  MatchEnvPtrVec deconv_envs;
  if (peak_list.size() > 0) {
    std::sort(peak_list.begin(), peak_list.end(), Peak::cmpPosInc);
    FlashDeconvSingleSp deconv(peak_list, kMs2TolerancePpm * 1e-6, max_charge,
                               max_mass, avg_ptr);
    PeakGroupVec peak_groups = deconv.deconv();
    // The spectrum's noise floor (DeconvData::initMinInte's estimate, reused
    // here so the reported envelopes are trimmed the way DeconvMs2Process
    // trims its own -- see toMatchEnvs); 0 when TopFD's noise estimation is
    // off, which also disables the trim, matching the standard pipeline.
    double min_inte = 0;
    double bound_inte = 0;
    if (topfd_para_ptr->isEstimateMinInte()) {
      bound_inte = baseline_util::getBaseLine(peak_list);
      double sn_ratio = topfd_para_ptr->getMsTwoSnRatio();
      min_inte = sn_ratio < 1.0 ? bound_inte * sn_ratio : bound_inte;
    }
    deconv_envs =
        toMatchEnvs(peak_groups, peak_list, env_para_ptr, min_inte, bound_inte);
  }
  // 3. Add precursor information
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  std::sort(sp_feat_ptr_vec.begin(), sp_feat_ptr_vec.end(),
            SpecFeature::cmpPrecInteDec);
  PrecursorPtrVec prec_ptr_vec;
  for (std::size_t i = 0; i < sp_feat_ptr_vec.size(); i++) {
    int prec_id = i;
    int feat_id = sp_feat_ptr_vec[i]->getFracFeatureId();
    double mono_mz = sp_feat_ptr_vec[i]->getPrecMonoMz();
    int charge = sp_feat_ptr_vec[i]->getPrecCharge();
    double inte = sp_feat_ptr_vec[i]->getPrecInte();
    PrecursorPtr prec_ptr =
        std::make_shared<Precursor>(prec_id, feat_id, mono_mz, charge, inte);
    prec_ptr_vec.push_back(prec_ptr);
  }
  header_ptr->setPrecPtrVec(prec_ptr_vec);

  // 4. Write to msalign file
  DeconvMsPtr deconv_ms_ptr =
      match_env_util::getDeconvMsPtr(header_ptr, deconv_envs);
  if (topfd_para_ptr->isOutputAvgMass()) {
    deconv_ms_ptr = avg_mass_util::getAvgMassMs(deconv_ms_ptr);
  }
  std::thread::id thread_id = std::this_thread::get_id();
  int writer_id = pool_ptr->getId(thread_id);
  ms2_writer_ptr_vec[writer_id]->writeMs(deconv_ms_ptr);

  // 5. Buffer the deconvoluted spectrum for the SQLite output
  if (sql_buffer_ptr != nullptr) {
    Ms2SqlRecord record{ms_ptr, std::move(deconv_envs)};
    sql_buffer_ptr->bytes += estimateRecordBytes(record);
    sql_buffer_ptr->thread_records[writer_id].push_back(std::move(record));
  }
}

std::function<void()> geneMsTwoTask(
    const MzmlMsPtr& ms_ptr, const SpecFeaturePtrVec& feat_ptr_vec,
    const TopfdParaPtr& topfd_para_ptr, const AveraginePtr& avg_ptr,
    const EnvParaPtr& env_para_ptr,
    const MsAlignWriterPtrVec& ms2_writer_ptr_vec,
    const SimpleThreadPoolPtr& pool_ptr,
    const Ms2SqlBufferPtr& sql_buffer_ptr) {
  return [ms_ptr, feat_ptr_vec, topfd_para_ptr, avg_ptr, env_para_ptr,
          ms2_writer_ptr_vec, pool_ptr, sql_buffer_ptr]() {
    deconvMsTwo(ms_ptr, feat_ptr_vec, topfd_para_ptr, avg_ptr, env_para_ptr,
                ms2_writer_ptr_vec, pool_ptr, sql_buffer_ptr);
  };
}

}  // namespace flashdeconv_ms2_process

FlashDeconvMs2Process::FlashDeconvMs2Process(
    const TopfdParaPtr& topfd_para_ptr, const std::string& output_filename_ext)
    : topfd_para_ptr_(topfd_para_ptr),
      output_filename_ext_(output_filename_ext) {}

void FlashDeconvMs2Process::readSpecFeature(
    const std::string& feat_file_name,
    std::map<int, SpecFeaturePtrVec>& feat_map) {
  SpecFeatureReaderPtr sp_feat_reader =
      std::make_shared<SpecFeatureReader>(feat_file_name);
  SpecFeaturePtrVec sp_feat_ptr_vec = sp_feat_reader->readAllFeatures();
  sp_feat_reader = nullptr;
  for (std::size_t i = 0; i < sp_feat_ptr_vec.size(); i++) {
    SpecFeaturePtr feat_ptr = sp_feat_ptr_vec[i];
    feat_map[feat_ptr->getSpecId()].push_back(feat_ptr);
  }
}

void FlashDeconvMs2Process::process() {
  MzmlMsGroupReaderPtr reader_ptr = std::make_shared<MzmlMsGroupReader>(
      topfd_para_ptr_->getMzmlFileName(), topfd_para_ptr_->getPrecWindowWidth(),
      topfd_para_ptr_->getActivation(), topfd_para_ptr_->getFracId(),
      topfd_para_ptr_->isFaims(), topfd_para_ptr_->getFaimsVoltage(),
      topfd_para_ptr_->isMissingLevelOne());

  MzmlMsGroupPtr ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  if (ms_group_ptr == nullptr) {
    LOG_ERROR("No spectrum to read in mzML file!");
    return;
  }
  // The averagine table (read-only, shared by the worker threads) and the
  // envelope parameters used to match the theoretical envelopes to the peaks
  flashdeconv_ms2_process::AveraginePtr avg_ptr =
      std::make_shared<flashdeconv_ms2_process::Averagine>(
          topfd_para_ptr_->getMaxMass());
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>(
      topfd_para_ptr_->getMzError(), topfd_para_ptr_->getMaxMissPeakNum(),
      topfd_para_ptr_->isRunFilterByMz(),
      topfd_para_ptr_->getMassGroupBoundary());
  // One SQLite writer shared across the worker threads; created only when
  // SQLite output is enabled.
  MzmlMsSqlWriterPtr sql_writer_ptr =
      topfd_para_ptr_->isGeneSql()
          ? std::make_shared<MzmlMsSqlWriter>(topfd_para_ptr_->getSqlDb())
          : nullptr;
  int thread_num = topfd_para_ptr_->getThreadNum();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(thread_num);
  flashdeconv_ms2_process::Ms2SqlBufferPtr sql_buffer_ptr =
      sql_writer_ptr != nullptr
          ? std::make_shared<flashdeconv_ms2_process::Ms2SqlBuffer>(thread_num)
          : nullptr;

  auto wait_until_pool_idle = [&]() {
    while (pool_ptr->getQueueSize() > 0 ||
           pool_ptr->getIdleThreadNum() < thread_num) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  };
  auto flush_sql_buffer = [&]() {
    for (auto& thread_buffer : sql_buffer_ptr->thread_records) {
      for (const auto& record : thread_buffer) {
        sql_writer_ptr->writeMs2(record.ms_ptr, record.deconv_envs);
      }
      thread_buffer.clear();
    }
    sql_buffer_ptr->bytes = 0;
  };

  std::string output_base_name = topfd_para_ptr_->getOutputBaseName();
  std::string ms2_msalign_name = output_base_name + "_" + output_filename_ext_;
  MsAlignWriterPtrVec ms2_writer_ptr_vec;
  for (int i = 0; i < thread_num; i++) {
    MsAlignWriterPtr ms2_ptr = std::make_shared<MsAlignWriter>(
        ms2_msalign_name + "_" + std::to_string(i));
    ms2_writer_ptr_vec.push_back(ms2_ptr);
  }
  std::map<int, SpecFeaturePtrVec> feat_map;
  SpecFeaturePtrVec missing_one_feat_list;
  if (topfd_para_ptr_->isMissingLevelOne()) {
    double max_mass = topfd_para_ptr_->getMaxMass();
    int max_charge = topfd_para_ptr_->getMaxCharge();
    double mono_mz = peak_util::compMz(max_mass, max_charge);
    SpecFeaturePtr feat_ptr =
        std::make_shared<SpecFeature>(mono_mz, max_charge);
    missing_one_feat_list.push_back(feat_ptr);
  } else {
    std::string feat_file_name = output_base_name + "_ms2.feature";
    readSpecFeature(feat_file_name, feat_map);
  }
  int spec_cnt = 0;
  int total_spec_num = topfd_para_ptr_->getMs2ScanNum();
  while (ms_group_ptr != nullptr) {
    MzmlMsPtrVec ms_ptr_vec = ms_group_ptr->getMsTwoPtrVec();
    for (std::size_t i = 0; i < ms_ptr_vec.size(); i++) {
      while (pool_ptr->getQueueSize() >=
             static_cast<std::size_t>(thread_num) * 2) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }
      spec_cnt++;
      MzmlMsPtr ms_ptr = ms_ptr_vec[i];
      SpecFeaturePtrVec sp_feat_ptr_vec;
      if (topfd_para_ptr_->isMissingLevelOne()) {
        sp_feat_ptr_vec = missing_one_feat_list;
      } else {
        auto feat_it = feat_map.find(ms_ptr->getMsHeaderPtr()->getSpecId());
        if (feat_it != feat_map.end()) {
          SpecFeaturePtrVec feat_list = feat_it->second;
          std::sort(feat_list.begin(), feat_list.end(),
                    SpecFeature::cmpPrecInteDec);
          sp_feat_ptr_vec.push_back(feat_list[0]);
          double first_inte = feat_list[0]->getPrecInte();
          for (std::size_t k = 1; k < feat_list.size(); k++) {
            if (feat_list[k]->getPrecInte() >=
                first_inte * topfd_para_ptr_->getPrecInteCutoffRatio()) {
              sp_feat_ptr_vec.push_back(feat_list[k]);
            }
          }
        }
      }
      pool_ptr->enqueue(flashdeconv_ms2_process::geneMsTwoTask(
          ms_ptr, sp_feat_ptr_vec, topfd_para_ptr_, avg_ptr, env_para_ptr,
          ms2_writer_ptr_vec, pool_ptr, sql_buffer_ptr));
      std::string msg = flashdeconv_ms2_process::updateMsTwoMsg(
          ms_ptr->getMsHeaderPtr(), spec_cnt, total_spec_num);
      std::cout << "\r" << msg << std::flush;
      if (sql_writer_ptr != nullptr &&
          sql_buffer_ptr->bytes.load() >=
              flashdeconv_ms2_process::SQL_BUFFER_FLUSH_BYTES) {
        wait_until_pool_idle();
        flush_sql_buffer();
      }
    }
    ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  }
  pool_ptr->shutDown();
  if (sql_writer_ptr != nullptr) {
    flush_sql_buffer();
    sql_writer_ptr->flush();
  }
  for (int i = 0; i < thread_num; i++) {
    ms2_writer_ptr_vec[i] = nullptr;
  }
  std::string para_str =
      topfd_para_ptr_->getParaStr("#", "\t", output_filename_ext_);
  MsalignThreadMergePtr ms2_merge_ptr = std::make_shared<MsalignThreadMerge>(
      output_filename_ext_, topfd_para_ptr_->getThreadNum(),
      output_filename_ext_, output_base_name, para_str);
  ms2_merge_ptr->process();
  std::string ms2_prefix = file_util::absoluteName(output_base_name) + "_" +
                           output_filename_ext_ + "_";
  file_util::cleanPrefix(output_base_name, ms2_prefix);
  std::cout << std::endl;
}

}  // namespace toppic
