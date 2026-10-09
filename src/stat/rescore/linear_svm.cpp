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

#include "stat/rescore/linear_svm.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <numeric>
#include <random>

namespace toppic {

namespace linear_svm {

double Model::score(const std::vector<double>& x) const {
  double s = bias;
  for (size_t j = 0; j < weights.size(); j++) {
    s += weights[j] * x[j];
  }
  return s;
}

Model train(const std::vector<std::vector<double>>& x,
            const std::vector<bool>& positive, double c, unsigned seed,
            int max_epoch, double eps) {
  size_t n = x.size();
  size_t d = n == 0 ? 0 : x[0].size();
  Model model;
  model.weights.assign(d, 0.0);
  if (n == 0) {
    return model;
  }
  size_t n_pos = std::count(positive.begin(), positive.end(), true);
  size_t n_neg = n - n_pos;
  // the dual's diagonal: x.x + 1 (the bias feature) + 1 / (2 C_i)
  std::vector<double> y(n);
  std::vector<double> diag(n);
  std::vector<double> q_bar(n);
  for (size_t i = 0; i < n; i++) {
    y[i] = positive[i] ? 1.0 : -1.0;
    size_t n_class = positive[i] ? n_pos : n_neg;
    double c_i = c * static_cast<double>(n) / (2.0 * n_class);
    diag[i] = 0.5 / c_i;
    double xx = 1.0;
    for (size_t j = 0; j < d; j++) {
      xx += x[i][j] * x[i][j];
    }
    q_bar[i] = xx + diag[i];
  }
  std::vector<double> alpha(n, 0.0);
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  // std::shuffle's algorithm is not fixed by the standard; a Fisher-Yates
  // shuffle on the generator's raw output gives the same order everywhere
  std::mt19937 rng(seed);
  for (int epoch = 0; epoch < max_epoch; epoch++) {
    for (size_t i = n - 1; i > 0; i--) {
      size_t k = rng() % (i + 1);
      std::swap(order[i], order[k]);
    }
    double pg_max = -std::numeric_limits<double>::infinity();
    double pg_min = std::numeric_limits<double>::infinity();
    for (size_t t = 0; t < n; t++) {
      size_t i = order[t];
      double g = y[i] * model.score(x[i]) - 1.0 + diag[i] * alpha[i];
      double pg = g;
      if (alpha[i] == 0.0) {
        pg = std::min(g, 0.0);
      }
      pg_max = std::max(pg_max, pg);
      pg_min = std::min(pg_min, pg);
      if (pg != 0.0) {
        double old_alpha = alpha[i];
        alpha[i] = std::max(alpha[i] - g / q_bar[i], 0.0);
        double delta = (alpha[i] - old_alpha) * y[i];
        for (size_t j = 0; j < d; j++) {
          model.weights[j] += delta * x[i][j];
        }
        model.bias += delta;
      }
    }
    if (pg_max - pg_min <= eps) {
      break;
    }
  }
  return model;
}

}  // namespace linear_svm

}  // namespace toppic
