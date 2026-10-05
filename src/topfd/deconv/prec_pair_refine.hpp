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

#ifndef TOPPIC_TOPFD_DECONV_PREC_PAIR_REFINE_HPP_
#define TOPPIC_TOPFD_DECONV_PREC_PAIR_REFINE_HPP_

#include "ms/env/match_env.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/precursor.hpp"

namespace toppic {

// Precursor masses from complementary fragment pairs (topfd's precursor
// isotope correction, on by default, --skip-ms2-prec-pair-refine): in HCD a b ion and its complementary y ion (b =
// residues, y = residues + water) sum to the proteoform's neutral mass
// whatever the proteoform's truncations and modifications, so the sums of
// pairs of deconvoluted fragment masses cluster at the precursor mass.
// Research switches through TOPFD_PREC_REFINE_PARA ("name=value,...").
namespace prec_pair_refine {

// Low-resolution MS1 (-R), which gives a precursor mass only to about 10
// ppm: for every precursor of a spectrum the pair sums within +-window Da
// of its mass are collected and the densest cluster of tol_ppm found; a
// precursor with at least min_pairs pairs in it takes the cluster's median
// as its mass. When the first precursor (the one TopPIC searches) has no
// such cluster and another has, that one becomes the first (the fragments
// are mostly its own).
void refine(const MatchEnvPtrVec& frag_envs, PrecursorPtrVec& precs);

// High-resolution MS1: the precursor mass is accurate to a few ppm but its
// monoisotopic peak can be wrong by whole isotopes. For every precursor
// (mass P) the pair sums of the fragment masses with an EnvCNN score of at
// least iso_min_envcnn are counted within iso_tol_ppm of P + j isotopes,
// j = -iso_range..iso_range. The precursor moves by the j with the most
// pairs (P + j isotopes, its m/z accuracy kept) when they number at least
// iso_min_pairs, more than at P and at least iso_ratio times the count at
// P, and, with ms1_peaks (the MS1 spectrum before the MS/MS spectrum, may
// be empty) and iso_ms1_check, when the averagine fit of the precursor in
// the MS1 spectrum at P beats the fit at P + j by less than
// iso_ms1_margin (cosines summed over the charges z-1..z+1). Precursors
// keep their order.
void correctIsotopes(const MatchEnvPtrVec& frag_envs,
                     const PeakPtrVec& ms1_peaks, PrecursorPtrVec& precs);

}  // namespace prec_pair_refine

}  // namespace toppic

#endif
