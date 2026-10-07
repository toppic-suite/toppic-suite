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

#ifndef TOPPIC_TOPFD_DECONV_PTCR_MS2_DECONV_HPP_
#define TOPPIC_TOPFD_DECONV_PTCR_MS2_DECONV_HPP_

#include <memory>
#include <string>
#include <vector>

#include "ms/mzml/mzml_ms.hpp"
#include "ms/spec/ms_header.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

// The deconvolution of PTCR spectra: MS/MS spectra of proton transfer charge
// reduction (PTCR) of the isolated ions without dissociation (activation
// PTCR, see PwMsReader). Proton transfer only removes charges, so every ion
// of such a spectrum is a product of an ion isolated in the precursor window
// [lo, hi] at some charge z0, and appears at z0, z0 - 1, ... with the mass
// and the isotope distribution of the isolated ion: one charge ladder per
// isolated species, as in an MS1 spectrum but bounded by the isolation.
//   1. The greedy MS1 search (greedy_ms1_deconv::deconvSeries; every charge
//      of a mass at one monoisotopic mass, chosen on the isotope intensities
//      summed over the charges) with the PTCR bound: charge z of mass M only
//      when some z0 >= z puts the envelope of M (its isotopes of at least
//      iso_min_rel of the apex) in the isolation window. The seed threshold
//      uses the MS/MS signal-to-noise ratio (-s).
//   2. One species per mass: its envelopes (EnvCNN >= min_envcnn) merged;
//      the mass refined from the matched isotope peaks of all its charges
//      (the intensity-weighted median of their mass errors), the summed
//      intensity, the charges.
//   3. The precursor: of the species with an envelope at the selected ion's
//      charge whose isotopes cover the selected ion's m/z, the most intense;
//      otherwise the most intense species that could have been isolated at
//      that charge; otherwise the most intense species; with none, the
//      selected ion itself.
// The research switches are in PtcrPara (cpp), set with the environment
// variable TOPFD_PTCR_PARA ("name=value,...").
namespace ptcr_ms2_deconv {

// One species of a PTCR spectrum: a mass with its charge ladder
struct Species {
  // the monoisotopic mass (refined from the matched peaks)
  double mono_mass = 0;
  // the mass of the search, before the refinement
  double search_mass = 0;
  // the summed intensity of its envelopes
  double inte = 0;
  // the best EnvCNN score of its envelopes
  double best_envcnn = 0;
  // the charges found, ascending, with the envelope intensity and EnvCNN
  // score of each
  std::vector<int> charges;
  std::vector<double> charge_intes;
  std::vector<double> charge_envcnns;
  // the charges z0 at which the mass's envelope overlaps the isolation
  // window (0: none, or no window)
  int iso_charge_min = 0;
  int iso_charge_max = 0;
  // the isotope peaks matched over all charges, and the weighted median
  // absolute deviation of their mass errors (ppm)
  int match_peak_num = 0;
  double mass_mad_ppm = 0;
  // the precursor of the spectrum
  bool is_prec = false;
};

// A deconvoluted PTCR spectrum: its header, with the precursor set, and its
// species, most intense first
struct PtcrSpec {
  MsHeaderPtr header_ptr;
  std::vector<Species> species;
};

using PtcrSpecPtr = std::shared_ptr<PtcrSpec>;
using PtcrSpecPtrVec = std::vector<PtcrSpecPtr>;

// Deconvolute one PTCR spectrum (MsHeader::isPtcr). Sets the precursor of
// the spectrum's header.
PtcrSpecPtr deconv(const MzmlMsPtr& ms_ptr, const TopfdParaPtr& topfd_para_ptr);

// Write the PTCR spectra, sorted by spectrum id: msalign_name in the msalign
// format (one line per species: mass, intensity, highest charge, best EnvCNN
// score; with -A average masses), tsv_name one line per species.
void write(PtcrSpecPtrVec spec_ptrs, const std::string& msalign_name,
           const std::string& tsv_name, const std::string& para_str,
           bool avg_mass);

}  // namespace ptcr_ms2_deconv

}  // namespace toppic

#endif
