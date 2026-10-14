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

#ifndef TOPPIC_STAT_RESCORE_LINEAR_SVM_HPP_
#define TOPPIC_STAT_RESCORE_LINEAR_SVM_HPP_

#include <vector>

namespace toppic {

namespace linear_svm {

// A linear support vector machine with an L2 penalty and the squared hinge
// loss (L2-loss SVM: LIBLINEAR's -s 1, scikit-learn's LinearSVC default),
// trained by dual coordinate descent (Hsieh et al., ICML 2008). The bias is
// one more feature of value 1, penalized with the weights, as LIBLINEAR and
// LinearSVC (intercept_scaling 1) do. The C of a sample is c x n / (2 x
// the size of its class) ("balanced" class weights).
struct Model {
  std::vector<double> weights;
  double bias = 0;

  double score(const std::vector<double>& x) const;
};

// x: one feature vector per sample; positive: the class of each sample.
// The samples are visited in an order shuffled with the seed, so the result
// is repeatable. Stops when the projected gradients of an epoch span at
// most eps, or after max_epoch epochs.
Model train(const std::vector<std::vector<double>>& x,
            const std::vector<bool>& positive, double c, unsigned seed,
            int max_epoch = 1000, double eps = 1e-3);

}  // namespace linear_svm

}  // namespace toppic

#endif
