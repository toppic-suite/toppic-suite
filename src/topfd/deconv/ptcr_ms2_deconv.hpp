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

#include "ms/env/match_env.hpp"
#include "ms/feature/spec_feature.hpp"
#include "ms/mzml/mzml_ms.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/spec/precursor.hpp"
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
//   3. The precursor: a species with an envelope at the selected ion's
//      charge whose monoisotopic m/z is the selected ion's (sel_ppm; Thermo
//      converters write the instrument's monoisotopic m/z as the selected
//      ion when it has one); otherwise of the species with an envelope at
//      that charge whose m/z range (from the monoisotopic peak to the last
//      isotope of at least iso_min_rel of the apex) covers the selected
//      ion's m/z, the most intense; otherwise the most intense species that
//      could have been isolated at that charge; otherwise the most intense
//      species; with none, the selected ion itself.
// With MS1 spectra, a species matching a feature of the spectrum's precursor
// window (10 ppm, also one isotope apart) gets the feature's id, and the
// precursor its id and intensity.
//
// Fragment PTCR spectra (MsHeader::isFragPtcr: a dissociation, then proton
// transfer charge reduction of the fragments, e.g. the MS3 scan "X@cid60
// Y@ptr10"): every fragment appears as a charge ladder from its charge z0,
// at which it lay in the isolation window of the proton transfer stage, and
// no fragment has more charges or mass than its precursor. They are
// deconvoluted as other MS/MS spectra (deconv_ms2_process::deconvMsTwo,
// --ms2-deconv) and written to the MS/MS msalign file with the
// dissociation's activation. The selected ion bounds the fragments
// (getSelectedIonPrec; with MS1 features, the larger of its bounds and the
// features'), and without MS1 features of the precursor window it is their
// precursor. The charge-series search above, with that bound,
// z0 <= the precursor's charge and masses below the precursor's
// (deconvFragments; switch frag=2), matched fewer theoretical fragments of
// AAV2 capsid proteins than the MS/MS search (toppic_test experiment 136,
// amendment 1).
//
// The dissociation of a charge-reduced precursor (MsHeader::isReducedPrec:
// proton transfer of the ions isolated at X, then the isolation of one
// charge-reduced ion at Y and its dissociation, e.g. the MS3 scan "X@ptr5
// Y@hcd80" after the PTCR MS2 scan of X) is an MS/MS spectrum of the
// dissociation's activation whose precursor is the species of the PTCR
// products at Y (getReducedPrec, from the parent PTCR spectrum; without it
// getReducedPrecFeatures).
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
  // the MS1 feature the species matches (-1: none) and its mass
  int feature_id = -1;
  double feature_mass = 0;
  // the envelopes of the species, one per charge
  MatchEnvPtrVec envs;
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
// the spectrum's header. feat_ptrs: the MS1 features of the spectrum's
// precursor window (none without MS1 spectra).
PtcrSpecPtr deconv(const MzmlMsPtr& ms_ptr, const TopfdParaPtr& topfd_para_ptr,
                   const SpecFeaturePtrVec& feat_ptrs);

// The precursor of the dissociation of a charge-reduced precursor
// (MsHeader::isReducedPrec), from the PTCR spectrum parent_ptr of the same
// isolation (the PTCR MS2 scan of the cycle): its species are deconvoluted
// again (the parent's header is read, not changed; feat_ptrs, the MS1
// features of the window, give them their feature ids), and the precursor
// is the most intense species with an envelope in the parent spectrum
// whose m/z range (monoisotopic peak to the last isotope of >= iso_min_rel)
// contains the target of the dissociation stage, else the most intense one
// whose envelope at a charge up to its isolation charge contains it: its
// mass at that charge, its intensity (its feature's when it has one), its
// feature id. nullptr without a parent or without such a species.
PrecursorPtr getReducedPrec(const MzmlMsPtr& parent_ptr,
                            const MsHeaderPtr& header_ptr,
                            const TopfdParaPtr& topfd_para_ptr,
                            const SpecFeaturePtrVec& feat_ptrs);

// The MS1 features of the window of the dissociation of a charge-reduced
// precursor (MsHeader::isReducedPrec) whose envelope at a charge up to
// theirs contains the target of the dissociation stage, in their order
SpecFeaturePtrVec getReducedPrecFeatures(const MsHeaderPtr& header_ptr,
                                         const SpecFeaturePtrVec& feat_ptrs);

// Deconvolute one fragment PTCR spectrum (MsHeader::isFragPtcr) with the
// charge-series search (switch frag=2; not the default). max_mass
// and max_charge bound the fragments (the precursors' largest mass and
// charge); feat_ptrs are the MS1 features of the precursor window (none
// without MS1 spectra). Sets the spectrum's precursors: the features';
// without them the precursor species of the spectrum (the surviving
// precursor, charge-reduced, at the selected ion's charge and m/z), else
// the selected ion. The species are the fragments, most intense first,
// without the masses within frag_prec_margin of the precursor mass or
// above it.
PtcrSpecPtr deconvFragments(const MzmlMsPtr& ms_ptr,
                            const TopfdParaPtr& topfd_para_ptr,
                            double max_mass, int max_charge,
                            const SpecFeaturePtrVec& feat_ptrs);

// The fragment species of a fragment PTCR spectrum as an MS/MS spectrum:
// one mass per species (intensity, highest charge, best EnvCNN score)
DeconvMsPtr getDeconvMs(const PtcrSpecPtr& spec_ptr);

// How fragment PTCR spectra are deconvoluted (research switch frag of
// TOPFD_PTCR_PARA): 1 (default) as other MS/MS spectra, the fragments
// bounded by the selected ion too, which is the precursor when the precursor
// window has no MS1 feature (getSelectedIonPrec); 0 as other MS/MS spectra
// (-o: the precursor of -c and -m); 2 deconvFragments
int getFragMode();

// The selected ion of an MS/MS spectrum (its m/z, else the isolation
// window's target, and its charge) as its precursor, nullptr without a
// charge. Sets max_mass and max_charge, the bounds of its fragments: the
// mass of the window's upper end at the selected charge, and that charge
// (-c and -m do not cap them).
PrecursorPtr getSelectedIonPrec(const MsHeaderPtr& header_ptr,
                                double& max_mass, int& max_charge);

// Write the species of PTCR or fragment PTCR spectra, one line per species,
// to tsv_name
void writeTsv(PtcrSpecPtrVec spec_ptrs, const std::string& tsv_name);

// Write the PTCR spectra, sorted by spectrum id: msalign_name in the msalign
// format (one line per species: mass, intensity, highest charge, best EnvCNN
// score; with -A average masses), tsv_name one line per species.
void write(PtcrSpecPtrVec spec_ptrs, const std::string& msalign_name,
           const std::string& tsv_name, const std::string& para_str,
           bool avg_mass);

}  // namespace ptcr_ms2_deconv

}  // namespace toppic

#endif
