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

#ifndef TOPPIC_STAT_RESCORE_RESCORE_MODEL_HPP_
#define TOPPIC_STAT_RESCORE_RESCORE_MODEL_HPP_

#include <vector>

namespace toppic {

namespace rescore_model {

struct Para {
  int fold_num = 3;
  int iter_num = 5;
  // the linear SVM's C
  double c = 0.1;
  // the FDR of the positives of every iteration
  double train_fdr = 0.01;
  unsigned seed = 1;
  // the ridge penalty of the retention time model
  double rt_alpha = 1.0;
};

// The input of one rescoring: one row per PrSM.
struct Data {
  std::vector<std::vector<double>> features;
  std::vector<bool> decoy;
  // the spectrum id of each PrSM: the PrSMs of a spectrum share a fold
  std::vector<int> spec_ids;
  // the column of features ranking the first iteration (higher is better)
  int init_col = 0;
  // the retention time model: a composition vector and the retention time
  // per PrSM; empty for none
  std::vector<std::vector<double>> comp;
  std::vector<double> rt;
};

// The targets passing the FDR (decoys with a score at least as good /
// targets at least as good, TopPIC's ratio) on a higher-is-better score;
// threshold: the smallest passing score, NaN when none passes.
std::vector<bool> passFdr(const std::vector<double>& scores,
                          const std::vector<bool>& decoy, double fdr,
                          double* threshold);

// Percolator's semi-supervised algorithm with cross-validation (Kall et al.,
// Nat. Methods 2007): the spectra are split into folds at random; for each
// fold a linear SVM is trained on the other folds, para.iter_num times,
// with the targets passing para.train_fdr by the current score as the
// positives and all decoys as the negatives (the first iteration ranks by
// the column init_col), and scores the fold; the fold scores are put on one
// scale (the training folds' threshold at para.train_fdr -> 0, their median
// decoy -> -1). With a retention time model, |rt - predicted rt| is one
// more feature, the prediction (ridge regression on comp) fitted in each
// fold on the training folds' targets passing para.train_fdr by init_col.
// Returns one score per row, higher is better, and false in ok when a fold
// had fewer than 10 positives (the scores are then those of init_col).
// mean_weights: the folds' mean weights on the standardized features (the
// retention time residual last).
std::vector<double> rescore(const Data& data, const Para& para, bool* ok,
                            std::vector<double>* mean_weights);

}  // namespace rescore_model

}  // namespace toppic

#endif
