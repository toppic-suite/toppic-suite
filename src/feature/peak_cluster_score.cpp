//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <fstream>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/util/file_util.hpp"
#include "feature/peak_cluster_score.hpp"

namespace toppic {

void loadTable(std::string file_name, std::vector<std::vector<double>> &table, 
               int row_num_, int col_num_) {
  std::ifstream input_;
  input_.open(file_name.c_str(), std::ios::in);
  if (!input_.is_open()) {
    LOG_ERROR("table file  " << file_name << " does not exist.");
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::vector<std::string> strs;
  for (int i = 0; i < row_num_; i++) {
    std::getline(input_, line);     
    strs = str_util::split(line, "\t");
    std::vector<double> one_row;
    for (int j = 0; j < col_num_; j++) {
      double v = std::stod(strs[j]);
      one_row.push_back(v);
    }
    table.push_back(one_row);
  }
}

PeakClusterScore::PeakClusterScore(std::string &dir, double thresh) {
  score_thresh_ = thresh;

  mass_bins_.resize(row_num_, 0.0);
  int idx = 0;
  for (int m = 800; m <= 30000; m += 1000) {
    mass_bins_[idx] = m;
    idx++;
  }
  std::string full_dir = dir + file_util::getFileSeparator();

  loadTable(full_dir + "DistScore.tsv", dist_score_table_, row_num_, col_num_);
  loadTable(full_dir + "CorrScore.tsv", corr_score_table_, row_num_, col_num_);
  loadTable(full_dir + "IntScore.tsv", inte_score_table_, row_num_, col_num_);

  loadTable(full_dir + "SummedDistScore.tsv", dist_sum_score_table_, row_num_, col_num_);
  loadTable(full_dir + "SummedCorrScore.tsv", corr_sum_score_table_, row_num_, col_num_);
  loadTable(full_dir + "SummedIntScore.tsv", inte_sum_score_table_, row_num_, col_num_);

  loadTable(full_dir + "AbuScore.tsv", inte_distr_score_table_, row_num_, col_num_);
  loadTable(full_dir + "XicCorrScore1.tsv", xic_score_table_1_, row_num_, col_num_);
  loadTable(full_dir + "XicCorrScore2.tsv", xic_score_table_2_, row_num_, col_num_);
}

double PeakClusterScore::getScore(PeakClusterPtr pc) {
  int mi = (int)std::round((pc->getMonoMass() - mass_bins_[0]) / (mass_bins_[1] - mass_bins_[0]));
  int max_idx = mass_bins_.size() - 1;
  mi = std::min(std::max(mi, 0), max_idx);
  
  double score = 0.0;

  for (int i = 0; i < 2; i++) {
    double inte_distr_value = pc->getInteDistr(i);
    int k = (int)std::min(std::max(std::round(inte_distr_value / 0.001), 0.0), col_num_ - 1.0);
    score += inte_distr_score_table_[mi][k];

    double dist_sum_value = std::min(pc->getSumDistScore(i), 1.0);
    k = (int)std::min(std::max(std::round(dist_sum_value / 0.001), 0.0), col_num_ - 1.0);
    score += dist_sum_score_table_[mi][k];

    double corr_sum_value = std::min(pc->getSumCorrScore(i), 1.0);
    k = (int)std::min(std::max(std::round(corr_sum_value / 0.001), 0.0), col_num_ - 1.0);
    score += corr_sum_score_table_[mi][k];

    double inte_sum_value = std::min(pc->getSumInteScore(i), 1.0);
    k = (int)std::min(std::max(std::round(inte_sum_value / 0.001), 0.0), col_num_ - 1.0);
    score += inte_sum_score_table_[mi][k];

    double best_dist_value = std::min(pc->getBestDistScore(i), 1.0);
    k = (int)std::min(std::max(std::round(best_dist_value / 0.001), 0.0), col_num_ - 1.0);
    score += dist_score_table_[mi][k];

    double best_corr_value = std::min(pc->getBestCorrScore(i), 1.0);
    k = (int)std::min(std::max(std::round(best_corr_value / 0.001), 0.0), col_num_ - 1.0);
    score += corr_score_table_[mi][k];

    double best_inte_value = std::min(pc->getBestInteScore(i), 1.0);
    k = (int)std::min(std::max(std::round(best_inte_value / 0.001), 0.0), col_num_ - 1.0);
    score += inte_score_table_[mi][k];

    double xic_value = std::min(pc->getXicCorrBetweenCharges(i), 1.0); 
    k = (int)std::min(std::max(std::round(xic_value / 0.001), 0.0), col_num_ - 1.0);
    double xic_score = (i == 0) ? xic_score_table_1_[mi][k] : xic_score_table_2_[mi][k];
    score += xic_score;
  }
  return score;
}

}
