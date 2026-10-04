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

#include "common/emass/emass.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <utility>

namespace toppic {

namespace emass {

namespace {

// Trailing peaks with an abundance of at most kLimit are pruned.
const double kLimit = 1e-100;

// The full convolution of u and v. Each output sums its products in the
// order of the longer series, as numpy.convolve does.
std::vector<double> convolveSeries(const std::vector<double> &u,
                                   const std::vector<double> &v) {
  const std::vector<double> &x = v.size() > u.size() ? v : u;
  const std::vector<double> &y = v.size() > u.size() ? u : v;
  std::vector<double> out(x.size() + y.size() - 1, 0.0);
  for (size_t i = 0; i < x.size(); i++) {
    for (size_t j = 0; j < y.size(); j++) {
      out[i + j] += x[i] * y[j];
    }
  }
  return out;
}

// Python's round(): halves to even, as the default rounding mode does.
double roundHalfEven(double x) { return std::nearbyint(x); }

// The sum of the terms with Neumaier's compensation, as Python's sum() of
// floats does from Python 3.12 (emass.py's monoisotopic masses).
double compensatedSum(const std::vector<double> &terms) {
  if (terms.empty()) {
    return 0.0;
  }
  double sum = terms[0];
  double c = 0.0;
  for (size_t i = 1; i < terms.size(); i++) {
    double x = terms[i];
    double t = sum + x;
    if (std::fabs(sum) >= std::fabs(x)) {
      c += (sum - t) + x;
    } else {
      c += (x - t) + sum;
    }
    sum = t;
  }
  if (c != 0 && std::isfinite(c)) {
    sum += c;
  }
  return sum;
}

}  // namespace

Formula averagineFormula(double mono_mass, AtomRounding rounding,
                         const Composition &composition) {
  std::vector<double> terms;
  for (int e = 0; e < kElementNum; e++) {
    terms.push_back(composition[e] *
                    getIsotopes(static_cast<Element>(e))[0].mass);
  }
  double n = mono_mass / compensatedSum(terms);
  Formula formula;
  for (int e = 0; e < kElementNum; e++) {
    formula[e] = composition[e] * n;
  }
  if (rounding == AtomRounding::kFractional) {
    return formula;
  }
  for (int e = 0; e < kElementNum; e++) {
    formula[e] = roundHalfEven(formula[e]);
  }
  if (rounding == AtomRounding::kRoundFillH) {
    std::vector<double> heavy;
    for (Element e : {kCarbon, kNitrogen, kOxygen, kSulfur}) {
      heavy.push_back(formula[e] * getIsotopes(e)[0].mass);
    }
    double rest = mono_mass - compensatedSum(heavy);
    formula[kHydrogen] = std::max(
        roundHalfEven(rest / getIsotopes(kHydrogen)[0].mass), 0.0);
  }
  return formula;
}

Formula averagineFormulaWithSulphur(double mono_mass, int sulphur_num,
                                    AtomRounding rounding,
                                    const Composition &composition) {
  double rest = mono_mass - sulphur_num * getIsotopes(kSulfur)[0].mass;
  Formula formula = averagineFormula(std::max(rest, 0.0), rounding,
                                     composition);
  formula[kSulfur] += sulphur_num;
  return formula;
}

double getMonoMass(const Formula &formula) {
  std::vector<double> terms;
  for (int e = 0; e < kElementNum; e++) {
    terms.push_back(formula[e] * getIsotopes(static_cast<Element>(e))[0].mass);
  }
  return compensatedSum(terms);
}

void getNeutronMoments(const Formula &formula, double *mean, double *var) {
  *mean = 0.0;
  *var = 0.0;
  for (int e = 0; e < kElementNum; e++) {
    const std::vector<Isotope> &iso = getIsotopes(static_cast<Element>(e));
    double sum = 0.0;
    double sum_k = 0.0;
    double sum_kk = 0.0;
    for (size_t k = 0; k < iso.size(); k++) {
      sum += iso[k].abundance;
      sum_k += k * iso[k].abundance;
      sum_kk += k * k * iso[k].abundance;
    }
    double mu = sum_k / sum;
    *mean += formula[e] * mu;
    *var += formula[e] * (sum_kk / sum - mu * mu);
  }
}

std::vector<double> seriesPower(const std::vector<double> &a, double n,
                                int k_max) {
  std::vector<double> b(k_max + 1, 0.0);
  b[0] = 1.0;
  int m = static_cast<int>(a.size()) - 1;
  for (int k = 1; k <= k_max; k++) {
    double s = 0.0;
    for (int j = 1; j <= std::min(k, m); j++) {
      s += ((n + 1) * j - k) * a[j] * b[k - j];
    }
    b[k] = s / k;
  }
  return b;
}

IsotopePattern::IsotopePattern(std::vector<double> a, std::vector<double> w)
    : a_(std::move(a)), w_(std::move(w)) {}

IsotopePattern IsotopePattern::getElementPattern(Element element) {
  const std::vector<Isotope> &iso = getIsotopes(element);
  double m0 = iso[0].mass;
  std::vector<double> a(iso.size());
  std::vector<double> w(iso.size());
  for (size_t k = 0; k < iso.size(); k++) {
    a[k] = iso[k].abundance;
    w[k] = iso[k].abundance > 0 ? iso[k].abundance * (iso[k].mass - m0) : 0.0;
  }
  return IsotopePattern(a, w);
}

IsotopePattern IsotopePattern::getElementPower(Element element, double n,
                                               int k_max) {
  const std::vector<Isotope> &iso = getIsotopes(element);
  double m0 = iso[0].mass;
  double p0 = iso[0].abundance;
  std::vector<double> a(iso.size());
  std::vector<double> w(iso.size());
  for (size_t k = 0; k < iso.size(); k++) {
    a[k] = iso[k].abundance / p0;
    w[k] = iso[k].abundance > 0
        ? iso[k].abundance * (iso[k].mass - m0) / p0 : 0.0;
  }
  std::vector<double> big = seriesPower(a, n, k_max);
  std::vector<double> w_n = convolveSeries(seriesPower(a, n - 1, k_max), w);
  w_n.resize(k_max + 1);
  for (double &x : w_n) {
    x *= n;
  }
  return IsotopePattern(big, w_n);
}

IsotopePattern IsotopePattern::convolve(const IsotopePattern &other) const {
  std::vector<double> a = convolveSeries(a_, other.a_);
  std::vector<double> w = convolveSeries(w_, other.a_);
  std::vector<double> w2 = convolveSeries(a_, other.w_);
  for (size_t k = 0; k < w.size(); k++) {
    w[k] += w2[k];
  }
  size_t end = a.size();
  while (end > 0 && !(a[end - 1] > kLimit)) {
    end--;
  }
  if (end == 0) {
    end = 1;
  }
  a.resize(end);
  w.resize(end);
  return IsotopePattern(a, w);
}

IsotopePattern IsotopePattern::power(long n) const {
  IsotopePattern result;
  IsotopePattern base = *this;
  while (n) {
    if (n & 1) {
      result = result.convolve(base);
    }
    n >>= 1;
    if (n) {
      base = base.convolve(base);
    }
  }
  return result;
}

void IsotopePattern::truncate(int k_max) {
  if (static_cast<int>(a_.size()) > k_max + 1) {
    a_.resize(k_max + 1);
    w_.resize(k_max + 1);
  }
}

std::vector<double> IsotopePattern::getCenters() const {
  std::vector<double> centers(a_.size());
  for (size_t k = 0; k < a_.size(); k++) {
    centers[k] = w_[k] / a_[k];
  }
  return centers;
}

IsotopeDistribution getDistribution(const Formula &formula,
                                    BelowOne below_one) {
  IsotopePattern pattern;
  int k_max = -1;
  for (int e = 0; e < kElementNum; e++) {
    Element element = static_cast<Element>(e);
    double n = formula[e];
    if (n == 0) {
      continue;
    }
    if (below_one == BelowOne::kMixture && n < 1) {
      IsotopePattern one = IsotopePattern::getElementPattern(element);
      std::vector<double> a = one.getAbundances();
      std::vector<double> w = one.getWeightedOffsets();
      for (size_t k = 0; k < a.size(); k++) {
        a[k] *= n;
        w[k] *= n;
      }
      a[0] += 1 - n;
      pattern = pattern.convolve(IsotopePattern(a, w));
    } else if (n == std::floor(n)) {
      pattern = pattern.convolve(
          IsotopePattern::getElementPattern(element).power(
              static_cast<long>(n)));
    } else {
      if (k_max < 0) {
        double mean;
        double var;
        getNeutronMoments(formula, &mean, &var);
        k_max = static_cast<int>(mean + 15 * std::sqrt(var) + 20);
      }
      pattern = pattern.convolve(
          IsotopePattern::getElementPower(element, n, k_max));
    }
  }
  if (k_max >= 0) {
    pattern.truncate(k_max);
  }
  IsotopeDistribution distr;
  distr.mono_mass = getMonoMass(formula);
  const std::vector<double> &a = pattern.getAbundances();
  double max_a = *std::max_element(a.begin(), a.end());
  distr.rel_abundances.resize(a.size());
  for (size_t k = 0; k < a.size(); k++) {
    distr.rel_abundances[k] = a[k] / max_a;
  }
  distr.mass_offsets = pattern.getCenters();
  return distr;
}

}  // namespace emass

}  // namespace toppic
