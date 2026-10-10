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

#ifndef TOPPIC_MS_CALIB_MASS_CALIB_MODEL_HPP_
#define TOPPIC_MS_CALIB_MASS_CALIB_MODEL_HPP_

#include <string>

namespace toppic {

// The m/z error model of one MS level, in ppm:
//   e = c0 + c_mz (m/z - mz0) / 1000 + c_rt_per_hour (rt - rt0) / 60
// with the retention time rt in minutes.
struct MassCalibLevel {
  double c0 = 0.0;
  double c_mz = 0.0;
  double c_rt_per_hour = 0.0;
};

// The two-pass m/z recalibration model (toppic_test experiment 152): one
// MassCalibLevel for the MS1 spectra (fitted on the precursor errors) and one
// for the MS/MS spectra (fitted on the fragment errors). A corrected m/z is
// m/z (1 - e 1e-6). The file format is a JSON object with the keys of the
// toppic_test scripts calib.py and recal_msalign.py (mz0, rt0_min, ms1, ms2,
// n_prec, n_frag, n_frag_prsms, robust_sd_prec, robust_sd_frag), so a model
// written by either can be read by the other.
class MassCalibModel {
 public:
  // the error (ppm) of an m/z value of an MS1 (ms_level 1) or MS/MS
  // (ms_level >= 2) spectrum at retention time rt_min (minutes)
  double getErrorPpm(int ms_level, double mz, double rt_min) const;

  // m/z (1 - e 1e-6)
  double correctMz(int ms_level, double mz, double rt_min) const;

  // the neutral mass of an ion of charge z whose m/z is corrected:
  // mass - e 1e-6 (mass + z proton), exact for mass = z (m/z - proton)
  double correctMass(int ms_level, double mass, int charge,
                     double rt_min) const;

  // writes the model as JSON
  void write(const std::string& file_name) const;

  // reads a model written by write() or by calib.py; exits with an error
  // message when a coefficient is missing
  static MassCalibModel read(const std::string& file_name);

  double mz0 = 800.0;
  double rt0_min = 0.0;
  MassCalibLevel ms1;
  MassCalibLevel ms2;
  // the fit's calibrant numbers and robust standard deviations (ppm), kept
  // for the report; not used by the correction
  int prec_num = 0;
  int frag_num = 0;
  int frag_prsm_num = 0;
  double prec_robust_sd = 0.0;
  double frag_robust_sd = 0.0;
};

}  // namespace toppic

#endif  // TOPPIC_MS_CALIB_MASS_CALIB_MODEL_HPP_
