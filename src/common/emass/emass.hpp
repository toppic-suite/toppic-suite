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

// Isotope distributions by the emass method (Rockwood and Haimi, JASMS
// 17:415, 2006; Rockwood, Van Orden and Smith, Anal Chem 67:2699, 1995),
// ported from toppic_test's emass.py (experiments 107, 108 and 111).
//
// A pattern is a list of peaks, one per number of extra neutrons, each with
// an abundance and the abundance-weighted mean mass of all isotopologues in
// it. Two patterns combine by convolution, the mass of each output peak being
// the abundance-weighted mean of the summed masses; the pattern of n atoms of
// an element is built by binary powers (squaring), and trailing peaks below
// 1e-100 are pruned. Masses are kept as offsets from the monoisotopic mass,
// so that 100 kDa molecules keep full precision.
//
// Fractional atom numbers: with the isotope polynomial P(x) = sum_k p_k x^k
// of an element and W(x) = sum_k p_k (m_k - m_0) x^k, n atoms have the
// abundances P(x)^n and the weighted mass offsets n P(x)^(n-1) W(x), which
// the convolutions compute for an integer n and which stay defined for any
// real n >= 0. The powers of a series are taken with J.C.P. Miller's
// recurrence (seriesPower).

#ifndef TOPPIC_COMMON_EMASS_EMASS_HPP_
#define TOPPIC_COMMON_EMASS_EMASS_HPP_

#include <array>
#include <vector>

#include "common/emass/isotope_data.hpp"

namespace toppic {

namespace emass {

// Numbers of C, H, N, O and S atoms, whole or fractional.
using Formula = std::array<double, kElementNum>;

// How averagineFormula makes the atom numbers from the scaled composition.
enum class AtomRounding {
  // the composition scaled to the mass, fractional atoms
  kFractional,
  // each of the five atom numbers rounded; the formula's mass is then up to
  // a few tens of Da from the mass asked for (TopFD's averagine table)
  kRoundEach,
  // C, N, O and S rounded and H filling the mass
  kRoundFillH,
};

// The composition scaled to the monoisotopic mass mono_mass.
Formula averagineFormula(double mono_mass,
                         AtomRounding rounding = AtomRounding::kRoundEach,
                         const Composition &composition = kHumanAveragine);

// The averagine formula of mono_mass with sulphur_num sulphur atoms: the
// composition, normally a sulphur-free one, scaled to mono_mass minus
// sulphur_num times the 32S mass (0 if that is negative) and rounded as
// averagineFormula does, plus sulphur_num S. With sulphur_num 0 it is
// averagineFormula.
Formula averagineFormulaWithSulphur(
    double mono_mass, int sulphur_num,
    AtomRounding rounding = AtomRounding::kRoundEach,
    const Composition &composition = kHumanNoSulphurAveragine);

// The monoisotopic mass of a formula.
double getMonoMass(const Formula &formula);

// Mean and variance of the number of extra neutrons of a formula.
void getNeutronMoments(const Formula &formula, double *mean, double *var);

// Coefficients 0..k_max of A(x)^n for a real n, A(x) = sum_j a[j] x^j with
// a[0] = 1 (J.C.P. Miller: k b_k = sum_j ((n + 1) j - k) a_j b_(k-j)); exact
// up to rounding for every k <= k_max, whatever is cut off above.
std::vector<double> seriesPower(const std::vector<double> &a, double n,
                                int k_max);

// Abundance a[k] and abundance x mass offset w[k] of peak k.
class IsotopePattern {
 public:
  // The pattern of nothing: one peak of abundance 1 at offset 0.
  IsotopePattern(): a_{1.0}, w_{0.0} {}

  IsotopePattern(std::vector<double> a, std::vector<double> w);

  // One atom of the element.
  static IsotopePattern getElementPattern(Element element);

  // n atoms of the element for a real n: peaks 0..k_max of P^n and
  // n P^(n-1) W, both divided by p_0^n (the patterns are used relative).
  static IsotopePattern getElementPower(Element element, double n, int k_max);

  // The pattern of the two together; trailing peaks with an abundance of at
  // most 1e-100 are removed.
  IsotopePattern convolve(const IsotopePattern &other) const;

  // n copies of the pattern, by binary powers.
  IsotopePattern power(long n) const;

  // Keeps peaks 0..k_max.
  void truncate(int k_max);

  const std::vector<double>& getAbundances() const { return a_; }

  const std::vector<double>& getWeightedOffsets() const { return w_; }

  // Mass offsets of the peak centers from the monoisotopic mass.
  std::vector<double> getCenters() const;

 private:
  std::vector<double> a_;
  std::vector<double> w_;
};

// How fewer than one atom of an element (0 < n < 1) is taken.
enum class BelowOne {
  // the formal series P(x)^n, whose tail has negative and unphysical terms
  kPower,
  // a fraction n of the molecules carrying one atom, (1 - n) + n P(x): the
  // same mean, equal to P(x)^n at n = 1, and a distribution at every n
  kMixture,
};

// The isotope distribution of a formula.
struct IsotopeDistribution {
  double mono_mass;
  // peak abundances relative to the highest peak
  std::vector<double> rel_abundances;
  // peak centers minus the monoisotopic mass
  std::vector<double> mass_offsets;
};

// Whole atom numbers use emass's binary powers; any other number
// getElementPower, up to 15 standard deviations above the mean number of
// extra neutrons (plus 20 peaks).
IsotopeDistribution getDistribution(const Formula &formula,
                                    BelowOne below_one = BelowOne::kPower);

}  // namespace emass

}  // namespace toppic

#endif
