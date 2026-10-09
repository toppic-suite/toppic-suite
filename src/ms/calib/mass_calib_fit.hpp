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

#ifndef TOPPIC_MS_CALIB_MASS_CALIB_FIT_HPP_
#define TOPPIC_MS_CALIB_MASS_CALIB_FIT_HPP_

#include <string>
#include <utility>
#include <vector>

#include "ms/calib/mass_calib_model.hpp"

namespace toppic {

// One calibrant: a precursor or a fragment mass of a confident PrSM, its m/z,
// the spectrum's retention time (minutes) and its error (ppm).
struct MassCalibrant {
  int scan = 0;
  double mass = 0.0;
  int charge = 0;
  double mz = 0.0;
  double rt_min = 0.0;
  double ppm = 0.0;
};

struct MassCalibrantSet {
  std::vector<MassCalibrant> prec;
  std::vector<MassCalibrant> frag;
  // the PrSMs of the table and those used (E-value, no shift, no variable PTM)
  int prsm_num = 0;
  int used_prsm_num = 0;
};

namespace mass_calib_fit {

// The precursor calibrant of a PrSM: the error (ppm of the proteoform mass)
// of the precursor mass against the proteoform mass, added to set->prec when
// the difference rounds to 0 isotopes of 1.00235 Da (and the charge is
// positive); returns whether it was added.
bool addPrecCalibrant(int scan, double rt_min, int charge, double prec_mass,
                      double form_mass, MassCalibrantSet* set);

// The fragment calibrants of a PrSM's spectrum: every mass (with its charge)
// within tol_ppm of a theoretical fragment mass (theo_masses sorted
// increasingly), with the error to the nearest one, added to set->frag;
// returns the number added.
int addFragCalibrants(int scan, double rt_min,
                      const std::vector<double>& theo_masses,
                      const std::vector<std::pair<double, int>>& masses,
                      double tol_ppm, MassCalibrantSet* set);

// The calibrants of a TopPIC PrSM table (a *_toppic_prsm_single.tsv file of a
// zero-shift search) and the MS/MS msalign file it searched: the PrSMs with
// an E-value <= max_evalue, no unexpected shift and no variable PTM.
// - Precursors: error = (precursor mass - proteoform mass) / proteoform mass
//   in ppm, for the PrSMs whose difference rounds to 0 isotopes (1.00235 Da);
//   the m/z is (mass + charge proton) / charge.
// - Fragments: the N- and C-terminal ion masses of the proteoform (residue
//   masses and fixed modifications from the base tables, the protein
//   N-terminal acetylation, the ion types of the PrSM's activation; b / y
//   for CID and HCD) matched to the spectrum's masses within tol_ppm, the
//   nearest one; each matched mass with its m/z and error.
// The proteoform string must hold only fixed modifications written as
// "(X)[Name]" on single residues and an optional "[Acetyl]-" prefix; other
// PrSMs are skipped. Exits with an error message on a missing column.
MassCalibrantSet readCalibrants(const std::string& prsm_tsv_file,
                                const std::string& ms2_msalign_file,
                                double max_evalue, double tol_ppm);

// The per-level model: MS1 from the precursors (c0 and the time term), MS/MS
// from the fragments (c0, the m/z term and the time term), each by robust
// least squares (10 rounds of Huber-type reweighting: weights
// sqrt(min(1, 1.5 s / |r|)), s = 1.4826 median |r|). rt0 is the median
// retention time of all calibrants; mz0 is 800.
MassCalibModel fitModel(const MassCalibrantSet& set);

// Robust least squares of y on the columns of x (row-major, cols columns);
// returns the coefficients and sets robust_sd (1.4826 median |residual|).
std::vector<double> robustFit(const std::vector<double>& x, int cols,
                              const std::vector<double>& y,
                              double* robust_sd);

// the median (the mean of the two middle values for an even count)
double median(std::vector<double> v);

}  // namespace mass_calib_fit

}  // namespace toppic

#endif  // TOPPIC_MS_CALIB_MASS_CALIB_FIT_HPP_
