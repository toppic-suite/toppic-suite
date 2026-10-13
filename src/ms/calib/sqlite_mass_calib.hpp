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

#ifndef TOPPIC_MS_CALIB_SQLITE_MASS_CALIB_HPP_
#define TOPPIC_MS_CALIB_SQLITE_MASS_CALIB_HPP_

#include <string>

#include "ms/calib/mass_calib_model.hpp"

namespace toppic {

namespace sqlite_mass_calib {

// Writes calib_db_name, TopFD's database ori_db_name with every m/z value and
// deconvoluted mass corrected with a calibration model (toppic --mass-calib,
// toppic_test experiment 154; calib_db_name is toppic's <base>_toppic.sqlite
// since experiment 156), for the search of the corrected spectra:
// - the MS1 model: ms1_peak.mz and ms1_env_peak.mz (m/z' = m/z (1 - e 1e-6)),
//   ms1_env.mono_mass and ms1_env.ref_mass (mass - e 1e-6 (mass + z proton),
//   e at the mass's m/z), and the ranges of the 3D levels of topfd -D (each
//   CONFIG row's MZMIN / MZMAX the range of its corrected level: since
//   toppic_test experiment 165 the rows of ms1_peak with LEVEL >= k, before
//   it the table PEAKS<k>, whose MZ is corrected too);
// - the MS/MS model: ms2_peak.mz, ms2_env_peak.mz, ms2_env.mono_mass and
//   ms2_env.ref_mass;
// each at its own spectrum's retention time (the spectrum table's
// retention_time in seconds, the PEAKS tables' RETENTIONTIME in ms). The
// isolation windows (ms2_spectrum.target_mz / begin_mz / end_mz), the
// intensities, retention times, ids and scores are copied unchanged, and
// toppic's own tables (prsm, proteoform, protein, prsm_mass_shift,
// prsm_protein_match, fasta_seq) are left out. A table mass_calib holds the
// model (MassCalibModel::toJson), the source database, the toppic version
// and the time. ori_db_name is only read; calib_db_name is replaced. The new
// database is written under a temporary name and renamed when it is
// complete, so an interrupted run leaves no partial file under its name.
// Exits with an error message when a database cannot be read or written.
void process(const MassCalibModel& model, const std::string& ori_db_name,
             const std::string& calib_db_name);

}  // namespace sqlite_mass_calib

}  // namespace toppic

#endif  // TOPPIC_MS_CALIB_SQLITE_MASS_CALIB_HPP_
