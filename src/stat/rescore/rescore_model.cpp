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

#include "stat/rescore/rescore_model.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <random>

#include "stat/rescore/linear_svm.hpp"

namespace toppic {

namespace rescore_model {

namespace {

using Matrix = std::vector<std::vector<double>>;

// rows of m selected by idx
Matrix selectRows(const Matrix& m, const std::vector<size_t>& idx) {
  Matrix out;
  out.reserve(idx.size());
  for (size_t i : idx) {
    out.push_back(m[i]);
  }
  return out;
}

template <typename T>
std::vector<T> selectValues(const std::vector<T>& v,
                            const std::vector<size_t>& idx) {
  std::vector<T> out;
  out.reserve(idx.size());
  for (size_t i : idx) {
    out.push_back(v[i]);
  }
  return out;
}

double median(std::vector<double> v) {
  if (v.empty()) {
    return 0.0;
  }
  size_t m = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + m, v.end());
  double hi = v[m];
  if (v.size() % 2 == 1) {
    return hi;
  }
  double lo = *std::max_element(v.begin(), v.begin() + m);
  return (lo + hi) / 2.0;
}

// solves a x = b by Gaussian elimination with partial pivoting (a small
// symmetric positive definite system)
std::vector<double> solve(Matrix a, std::vector<double> b) {
  size_t n = b.size();
  for (size_t col = 0; col < n; col++) {
    size_t piv = col;
    for (size_t r = col + 1; r < n; r++) {
      if (std::abs(a[r][col]) > std::abs(a[piv][col])) {
        piv = r;
      }
    }
    std::swap(a[col], a[piv]);
    std::swap(b[col], b[piv]);
    double p = a[col][col];
    if (p == 0.0) {
      continue;
    }
    for (size_t r = col + 1; r < n; r++) {
      double f = a[r][col] / p;
      if (f == 0.0) {
        continue;
      }
      for (size_t k = col; k < n; k++) {
        a[r][k] -= f * a[col][k];
      }
      b[r] -= f * b[col];
    }
  }
  std::vector<double> x(n, 0.0);
  for (size_t r = n; r-- > 0;) {
    double s = b[r];
    for (size_t k = r + 1; k < n; k++) {
      s -= a[r][k] * x[k];
    }
    x[r] = a[r][r] == 0.0 ? 0.0 : s / a[r][r];
  }
  return x;
}

// ridge regression with an unpenalized intercept: the weights and, last,
// the intercept
std::vector<double> ridgeFit(const Matrix& x, const std::vector<double>& y,
                             double alpha) {
  size_t n = x.size();
  size_t d = n == 0 ? 0 : x[0].size();
  std::vector<double> mean_x(d, 0.0);
  double mean_y = 0.0;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < d; j++) {
      mean_x[j] += x[i][j];
    }
    mean_y += y[i];
  }
  for (size_t j = 0; j < d; j++) {
    mean_x[j] /= std::max<size_t>(n, 1);
  }
  mean_y /= std::max<size_t>(n, 1);
  Matrix a(d, std::vector<double>(d, 0.0));
  std::vector<double> b(d, 0.0);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < d; j++) {
      double xj = x[i][j] - mean_x[j];
      b[j] += xj * (y[i] - mean_y);
      for (size_t k = j; k < d; k++) {
        a[j][k] += xj * (x[i][k] - mean_x[k]);
      }
    }
  }
  for (size_t j = 0; j < d; j++) {
    for (size_t k = 0; k < j; k++) {
      a[j][k] = a[k][j];
    }
    a[j][j] += alpha;
  }
  std::vector<double> w = solve(a, b);
  double intercept = mean_y;
  for (size_t j = 0; j < d; j++) {
    intercept -= w[j] * mean_x[j];
  }
  w.push_back(intercept);
  return w;
}

double ridgePredict(const std::vector<double>& w,
                    const std::vector<double>& x) {
  double s = w.back();
  for (size_t j = 0; j < x.size(); j++) {
    s += w[j] * x[j];
  }
  return s;
}

}  // namespace

std::vector<bool> passFdr(const std::vector<double>& scores,
                          const std::vector<bool>& decoy, double fdr,
                          double* threshold) {
  size_t n = scores.size();
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&scores](size_t a, size_t b) {
    return scores[a] > scores[b];
  });
  std::vector<bool> pass(n, false);
  *threshold = std::numeric_limits<double>::quiet_NaN();
  int n_decoy = 0;
  int n_target = 0;
  size_t i = 0;
  while (i < n) {
    // a group of equal scores: its decoys count against its targets
    size_t j = i;
    while (j < n && scores[order[j]] == scores[order[i]]) {
      if (decoy[order[j]]) {
        n_decoy++;
      }
      j++;
    }
    for (size_t k = i; k < j; k++) {
      if (!decoy[order[k]]) {
        n_target++;
        double ratio = std::min(1.0, static_cast<double>(n_decoy) / n_target);
        if (ratio <= fdr) {
          pass[order[k]] = true;
          *threshold = scores[order[k]];
        }
      }
    }
    i = j;
  }
  return pass;
}

std::vector<double> rescore(const Data& data, const Para& para, bool* ok,
                            std::vector<double>* mean_weights) {
  size_t n = data.features.size();
  std::vector<double> init(n);
  for (size_t i = 0; i < n; i++) {
    init[i] = data.features[i][data.init_col];
  }
  *ok = false;
  mean_weights->clear();
  if (n == 0) {
    return init;
  }
  bool use_rt = !data.comp.empty();
  // folds by spectrum, in increasing spectrum id order
  std::map<int, int> fold_of;
  for (int id : data.spec_ids) {
    fold_of[id] = 0;
  }
  std::mt19937 rng(para.seed);
  for (auto& f : fold_of) {
    f.second = static_cast<int>(rng() % para.fold_num);
  }
  std::vector<double> result(n, 0.0);
  size_t d = data.features[0].size() + (use_rt ? 1 : 0);
  std::vector<double> weight_sum(d, 0.0);
  for (int fold = 0; fold < para.fold_num; fold++) {
    std::vector<size_t> tr;
    std::vector<size_t> te;
    for (size_t i = 0; i < n; i++) {
      if (fold_of[data.spec_ids[i]] == fold) {
        te.push_back(i);
      } else {
        tr.push_back(i);
      }
    }
    if (te.empty()) {
      continue;
    }
    Matrix x_tr = selectRows(data.features, tr);
    Matrix x_te = selectRows(data.features, te);
    std::vector<bool> decoy_tr = selectValues(data.decoy, tr);
    std::vector<double> init_tr = selectValues(init, tr);
    double thr = 0;
    if (use_rt) {
      // the retention time model on the training folds' confident targets
      std::vector<bool> pos = passFdr(init_tr, decoy_tr, para.train_fdr, &thr);
      Matrix comp_pos;
      std::vector<double> rt_pos;
      for (size_t k = 0; k < tr.size(); k++) {
        if (pos[k]) {
          comp_pos.push_back(data.comp[tr[k]]);
          rt_pos.push_back(data.rt[tr[k]]);
        }
      }
      if (comp_pos.size() < 10) {
        return init;
      }
      std::vector<double> w = ridgeFit(comp_pos, rt_pos, para.rt_alpha);
      for (size_t k = 0; k < tr.size(); k++) {
        x_tr[k].push_back(
            std::abs(data.rt[tr[k]] - ridgePredict(w, data.comp[tr[k]])));
      }
      for (size_t k = 0; k < te.size(); k++) {
        x_te[k].push_back(
            std::abs(data.rt[te[k]] - ridgePredict(w, data.comp[te[k]])));
      }
    }
    // standardization on the training folds
    std::vector<double> mean(d, 0.0);
    std::vector<double> sd(d, 0.0);
    for (const auto& row : x_tr) {
      for (size_t j = 0; j < d; j++) {
        mean[j] += row[j];
      }
    }
    for (size_t j = 0; j < d; j++) {
      mean[j] /= x_tr.size();
    }
    for (const auto& row : x_tr) {
      for (size_t j = 0; j < d; j++) {
        sd[j] += (row[j] - mean[j]) * (row[j] - mean[j]);
      }
    }
    for (size_t j = 0; j < d; j++) {
      sd[j] = std::sqrt(sd[j] / x_tr.size());
      if (sd[j] == 0.0) {
        sd[j] = 1.0;
      }
    }
    auto standardize = [&mean, &sd, d](Matrix* m) {
      for (auto& row : *m) {
        for (size_t j = 0; j < d; j++) {
          row[j] = (row[j] - mean[j]) / sd[j];
        }
      }
    };
    standardize(&x_tr);
    standardize(&x_te);
    std::vector<double> score_tr = init_tr;
    linear_svm::Model model;
    bool trained = false;
    for (int it = 0; it < para.iter_num; it++) {
      std::vector<bool> pos = passFdr(score_tr, decoy_tr, para.train_fdr, &thr);
      Matrix x_fit;
      std::vector<bool> y_fit;
      size_t n_pos = 0;
      for (size_t k = 0; k < tr.size(); k++) {
        if (pos[k] || decoy_tr[k]) {
          x_fit.push_back(x_tr[k]);
          y_fit.push_back(pos[k]);
          if (pos[k]) {
            n_pos++;
          }
        }
      }
      if (n_pos < 10) {
        break;
      }
      model = linear_svm::train(x_fit, y_fit, para.c,
                                para.seed + 1000u * fold + it);
      trained = true;
      for (size_t k = 0; k < tr.size(); k++) {
        score_tr[k] = model.score(x_tr[k]);
      }
    }
    if (!trained) {
      return init;
    }
    // one scale for the folds
    std::vector<double> decoy_scores;
    for (size_t k = 0; k < tr.size(); k++) {
      if (decoy_tr[k]) {
        decoy_scores.push_back(score_tr[k]);
      }
    }
    passFdr(score_tr, decoy_tr, para.train_fdr, &thr);
    if (std::isnan(thr)) {
      thr = *std::max_element(score_tr.begin(), score_tr.end());
    }
    double med_d = median(decoy_scores);
    double scale = thr > med_d ? thr - med_d : 1.0;
    for (size_t k = 0; k < te.size(); k++) {
      result[te[k]] = (model.score(x_te[k]) - thr) / scale;
    }
    for (size_t j = 0; j < d; j++) {
      weight_sum[j] += model.weights[j] / para.fold_num;
    }
  }
  *ok = true;
  *mean_weights = weight_sum;
  return result;
}

}  // namespace rescore_model

}  // namespace toppic
