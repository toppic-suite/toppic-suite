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

#include "search/graphalign/graph_align.hpp"

#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "ms/factory/extend_ms_factory.hpp"
#include "search/diag/diag_header_util.hpp"
#include "search/diag/diagonal_util.hpp"
#include "search/graph/graph.hpp"
#include "search/graphalign/graph_align_processor.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

std::vector<int> getMinMaxProtDist(const DistList2D& dist_lists) {
  std::vector<int> min_dist, max_dist;
  // the range of all the distances, including those a cover left out
  for (size_t i = 0; i < dist_lists.size(); i++) {
    if (dist_lists[i].all_num_ == 0) continue;
    min_dist.push_back(dist_lists[i].all_min_);
    max_dist.push_back(dist_lists[i].all_max_);
  }
  std::sort(min_dist.begin(), min_dist.end());
  std::sort(max_dist.begin(), max_dist.end());
  std::vector<int> res;
  res.push_back(min_dist[0]);
  res.push_back(max_dist[max_dist.size() - 1]);
  return res;
}

// A score read into an int. The scores are integers, or -DBL_MAX for an entry
// the DP has not reached; the conversion of -DBL_MAX to int is undefined, and
// x86 gives INT_MIN, which is what the DP compares here. Every score below -1
// compares the same with the -1 the searches start from.
inline int toIntScore(double score) {
  if (score < static_cast<double>(std::numeric_limits<int>::min())) {
    return std::numeric_limits<int>::min();
  }
  return static_cast<int>(score);
}

GraphAlign::GraphAlign(const GraphAlignMngPtr& mng_ptr,
                       const ProteoGraphPtr& proteo_graph_ptr,
                       const SpecGraphPtr& spec_graph_ptr,
                       const FastaSeqPtr& fasta_seq_ptr) {
  LOG_DEBUG("Graph constructor start");
  mng_ptr_ = mng_ptr;
  proteo_graph_ptr_ = proteo_graph_ptr;
  spec_graph_ptr_ = spec_graph_ptr;
  fasta_seq_ptr_ = fasta_seq_ptr;

  // The distance lists of both graphs are in increasing order of distance,
  // one entry per distance (dist.hpp), so they are read in place: neither is
  // copied or sorted again.
  std::vector<int> cutoff =
      getMinMaxProtDist(proteo_graph_ptr_->getDistList2D());
  cutoff[0] -= mng_ptr_->getIntTolerance();
  cutoff[1] += mng_ptr_->getIntTolerance();

  const DistList& spec_dist = spec_graph_ptr_->getDistList();
  for (size_t i = 0; i < spec_dist.size(); i++) {
    if (spec_dist.dists_[i] >= cutoff[0] && spec_dist.dists_[i] <= cutoff[1]) {
      spec_dist_idxes_.push_back(static_cast<int>(i));
    }
  }

  pg_ = proteo_graph_ptr_->getMassGraphPtr();
  sg_ = spec_graph_ptr_->getMassGraphPtr();
  proteo_ver_num_ = num_vertices(*pg_.get());
  spec_ver_num_ = num_vertices(*sg_.get());
  shift_layer_num_ = mng_ptr_->n_unknown_shift_ + 1;
  mod_layer_num_ = mng_ptr_->max_known_mods_ + 1;

  // As toppic's Proteoform::setTruncNTermModPtr allows a truncation
  // N-terminal modification: on a proteoform starting at protein residue 3
  // or later (0-based 2; the first two are the protein N-terminal forms')
  // whose first amino acid is one the modification allows.
  ProteoformPtr db_form_ptr = proteo_graph_ptr_->getProteoformPtr();
  int seq_start = db_form_ptr->getStartPos();
  ResSeqPtr res_seq_ptr = db_form_ptr->getResSeqPtr();
  n_term_starts_.assign(mng_ptr_->n_term_int_masses_.size(),
                        std::vector<char>(proteo_ver_num_, 0));
  for (size_t t = 0; t < n_term_starts_.size(); t++) {
    const std::vector<AminoAcidPtr>& acid_ptrs = mng_ptr_->n_term_acid_ptrs_[t];
    for (int i = std::max(0, 2 - seq_start); i < res_seq_ptr->getLen(); i++) {
      AminoAcidPtr acid_ptr = res_seq_ptr->getResiduePtr(i)->getAminoAcidPtr();
      if (std::find(acid_ptrs.begin(), acid_ptrs.end(), acid_ptr) !=
          acid_ptrs.end()) {
        n_term_starts_[t][i] = 1;
      }
    }
  }
  LOG_DEBUG("Graph constructor end");
}

void GraphAlign::getConsistentPairs() {
  LOG_DEBUG("consistent pair start");
  int tole = mng_ptr_->getIntTolerance();
  LOG_DEBUG("Integer error tolerance " << tole);
  const DistList2D& prot_lists = proteo_graph_ptr_->getDistList2D();
  const DistList& spec_list = spec_graph_ptr_->getDistList();
  int min_dist = mng_ptr_->getIntMinConsistentDist();

  // Every pair of a proteoform distance and a spectrum distance within the
  // tolerance makes the pairs of their vertex pairs consistent. The search
  // runs twice, first counting the pairs of each cell and modification
  // number, then storing them, so that they are stored in the order they
  // are found.
  size_t key_num = static_cast<size_t>(proteo_ver_num_) * spec_ver_num_ *
                   mod_layer_num_;
  std::vector<int> cursors(key_num + 1, 0);
  for (int pass = 0; pass < 2; pass++) {
    for (size_t m = 0; m < prot_lists.size(); m++) {
      const DistList& prot_list = prot_lists[m];
      if (prot_list.size() == 0) continue;

      size_t prot_idx_min = 0, prot_idx_max = 0;
      for (size_t spec_idx = 0; spec_idx < spec_dist_idxes_.size();
           spec_idx++) {
        int sd = spec_dist_idxes_[spec_idx];
        int sp_dist = spec_list.dists_[sd];
        if (sp_dist < min_dist) continue;

        while (prot_idx_min < prot_list.size() &&
               prot_list.dists_[prot_idx_min] < sp_dist - tole) {
          prot_idx_min++;
        }

        if (prot_idx_min >= prot_list.size()) continue;

        prot_idx_max = std::max(prot_idx_min, prot_idx_max);

        while (prot_idx_max < prot_list.size() &&
               prot_list.dists_[prot_idx_max] <= sp_dist + tole) {
          prot_idx_max++;
        }

        for (size_t t = prot_idx_min; t < prot_idx_max; t++) {
          if (std::abs(sp_dist - prot_list.dists_[t]) > tole) continue;
          for (int k = prot_list.begins_[t]; k < prot_list.begins_[t + 1];
               k++) {
            const std::pair<int, int>& pg_pair = prot_list.pairs_[k];
            for (int sp = spec_list.begins_[sd]; sp < spec_list.begins_[sd + 1];
                 sp++) {
              const std::pair<int, int>& sp_pair = spec_list.pairs_[sp];
              size_t key =
                  getCellIdx(pg_pair.second, sp_pair.second) * mod_layer_num_ +
                  m;
              if (pass == 0) {
                cursors[key + 1]++;
              } else {
                cons_prev_cells_[cursors[key]++] = static_cast<int>(
                    getCellIdx(pg_pair.first, sp_pair.first));
              }
            }
          }
        }
      }
    }
    // Proteoform N-terminal variable PTMs: a distance of the spectrum's zero
    // vertex (the N-terminus) also matches a proteoform distance m shorter
    // from a vertex where the modification of mass m is allowed. The
    // alignment then starts at that vertex with the modification, which
    // takes one variable PTM layer more, as the protein N-terminal forms'
    // edges do; GraphPostProcessor reads it back from the variable shift
    // that begins at the first residue. The pair stores its cell c as -1 - c,
    // so that the DP knows the path has the modification.
    for (size_t t = 0; t < n_term_starts_.size(); t++) {
      int n_term_mass = mng_ptr_->n_term_int_masses_[t];
      const std::vector<char>& starts = n_term_starts_[t];
      for (size_t m = 0; m + 1 < prot_lists.size(); m++) {
        const DistList& prot_list = prot_lists[m];
        if (prot_list.size() == 0) continue;
        for (int sp = 1; sp < spec_ver_num_; sp++) {
          int dist = spec_graph_ptr_->getPeakDist(0, sp) - n_term_mass;
          if (dist < min_dist) continue;
          size_t t_idx = std::lower_bound(prot_list.dists_.begin(),
                                          prot_list.dists_.end(), dist - tole) -
                         prot_list.dists_.begin();
          for (; t_idx < prot_list.size() &&
                 prot_list.dists_[t_idx] <= dist + tole;
               t_idx++) {
            for (int k = prot_list.begins_[t_idx];
                 k < prot_list.begins_[t_idx + 1]; k++) {
              const std::pair<int, int>& pg_pair = prot_list.pairs_[k];
              if (!starts[pg_pair.first]) continue;
              size_t key =
                  getCellIdx(pg_pair.second, sp) * mod_layer_num_ + m + 1;
              if (pass == 0) {
                cursors[key + 1]++;
              } else {
                cons_prev_cells_[cursors[key]++] =
                    -1 - static_cast<int>(getCellIdx(pg_pair.first, 0));
              }
            }
          }
        }
      }
    }
    if (pass == 0) {
      for (size_t key = 0; key < key_num; key++) {
        cursors[key + 1] += cursors[key];
      }
      cons_begins_ = cursors;
      cons_prev_cells_.resize(cursors[key_num]);
    }
  }
  LOG_DEBUG("consistent pair end");
}

void GraphAlign::initTable() {
  LOG_DEBUG("init table start");
  double node_score = 1.0;
  size_t entry_num = static_cast<size_t>(proteo_ver_num_) * spec_ver_num_ *
                     shift_layer_num_ * mod_layer_num_;
  prev_edge_types_.assign(entry_num, GRAPH_ALIGN_TYPE_NULL);
  prev_edge_mod_nums_.assign(entry_num, 0);
  prev_nodes_.assign(entry_num, -1);
  best_scores_.assign(entry_num, -std::numeric_limits<double>::max());
  best_shift_nodes_.assign(entry_num, -1);
  best_shift_scores_.assign(entry_num, -std::numeric_limits<double>::max());
  trunc_paths_.assign(entry_num, 0);
  LOG_DEBUG("init table step 1");
  int init_aa_num = proteo_ver_num_;
  if (mng_ptr_->whole_protein_only_) {
    // the protein N-terminus, or its second residue when NME is allowed for
    // the protein (the acetylated forms start at the N-terminus)
    init_aa_num = std::min(proteo_graph_ptr_->isNme() ? 2 : 1, proteo_ver_num_);
  }
  for (int i = 0; i < init_aa_num; i++) {
    size_t cell = getCellIdx(i, 0);
    size_t entry = getEntryIdx(cell, 0, 0);
    // the start of an alignment: its score is the node score as an int
    prev_edge_types_[entry] = GRAPH_ALIGN_TYPE_NULL;
    prev_edge_mod_nums_[entry] = 0;
    prev_nodes_[entry] = -1;
    best_scores_[entry] = static_cast<int>(node_score);
    best_shift_scores_[entry] = 0;
    best_shift_nodes_[entry] = static_cast<int>(cell);
  }
  LOG_DEBUG("init table end");
}

// On equal scores, the DP keeps the path without a truncation N-terminal
// modification, so that a protein N-terminal form wins a tie against it; the
// paths are otherwise chosen as without the modifications.
int GraphAlign::compBestVariableNode(size_t cell, int s, int m,
                                     int& best_edge_mod_num,
                                     bool& best_trunc) {
  int best_prev_score = -1;
  best_edge_mod_num = -1;
  best_trunc = false;
  int best_prev_node = -1;
  for (int p = 0; p <= m; p++) {
    size_t key = cell * mod_layer_num_ + p;
    for (int q = cons_begins_[key]; q < cons_begins_[key + 1]; q++) {
      int prev_cell = cons_prev_cells_[q];
      bool trunc = prev_cell < 0;
      if (trunc) {
        prev_cell = -1 - prev_cell;
      }
      size_t prev_entry = getEntryIdx(prev_cell, s, m - p);
      int score = toIntScore(best_scores_[prev_entry]);
      trunc = trunc || trunc_paths_[prev_entry];
      if (score > best_prev_score ||
          (score == best_prev_score && best_trunc && !trunc)) {
        best_prev_score = score;
        best_prev_node = prev_cell;
        best_edge_mod_num = p;
        best_trunc = trunc;
      }
    }
  }
  return best_prev_node;
}

int GraphAlign::compBestShiftNode(int i, int j, int s, int m) {
  if (s == 0) {
    return -1;
  }
  int best_prev_score = -1;
  int best_prev_node = -1;
  size_t up_entry = getEntryIdx(getCellIdx(i - 1, j), s - 1, m);
  int score = toIntScore(best_shift_scores_[up_entry]);
  if (score > best_prev_score) {
    best_prev_score = score;
    best_prev_node = best_shift_nodes_[up_entry];
  }

  size_t left_entry = getEntryIdx(getCellIdx(i, j - 1), s - 1, m);
  score = toIntScore(best_shift_scores_[left_entry]);
  int left_node = best_shift_nodes_[left_entry];
  if (score > best_prev_score ||
      (score == best_prev_score && isTruncPath(best_prev_node, s - 1, m) &&
       !isTruncPath(left_node, s - 1, m))) {
    best_prev_score = score;
    best_prev_node = left_node;
  }
  return best_prev_node;
}

void GraphAlign::updateBestShiftNode(int i, int j, int s, int m) {
  size_t cell = getCellIdx(i, j);
  size_t entry = getEntryIdx(cell, s, m);
  size_t up_entry = getEntryIdx(getCellIdx(i - 1, j), s, m);
  size_t left_entry = getEntryIdx(getCellIdx(i, j - 1), s, m);
  // update best node
  double up_score = best_shift_scores_[up_entry];
  double left_score = best_shift_scores_[left_entry];
  if (up_score > left_score ||
      (up_score == left_score &&
       isTruncPath(best_shift_nodes_[left_entry], s, m) &&
       !isTruncPath(best_shift_nodes_[up_entry], s, m))) {
    best_shift_scores_[entry] = up_score;
    best_shift_nodes_[entry] = best_shift_nodes_[up_entry];
  } else {
    best_shift_scores_[entry] = left_score;
    best_shift_nodes_[entry] = best_shift_nodes_[left_entry];
  }

  if (best_scores_[entry] > best_shift_scores_[entry] ||
      (best_scores_[entry] == best_shift_scores_[entry] &&
       isTruncPath(best_shift_nodes_[entry], s, m) && !trunc_paths_[entry])) {
    best_shift_scores_[entry] = best_scores_[entry];
    best_shift_nodes_[entry] = static_cast<int>(cell);
  }
}

void GraphAlign::dp() {
  LOG_DEBUG("dp start");
  double new_score = 1.0;
  for (int i = 1; i < proteo_ver_num_; i++) {
    for (int j = 1; j < spec_ver_num_; j++) {
      size_t cell = getCellIdx(i, j);
      // compute for zero shift
      for (int s = 0; s <= mng_ptr_->n_unknown_shift_; s++) {
        for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
          int edge_mod_num;
          bool var_trunc;
          int best_var_node =
              compBestVariableNode(cell, s, m, edge_mod_num, var_trunc);
          double var_score;
          if (best_var_node < 0) {
            var_score = -std::numeric_limits<double>::max();
          } else {
            var_score =
                best_scores_[getEntryIdx(best_var_node, s, m - edge_mod_num)];
          }

          int best_shift_node = compBestShiftNode(i, j, s, m);
          double shift_score;
          bool shift_trunc = isTruncPath(best_shift_node, s - 1, m);
          if (best_shift_node >= 0) {
            shift_score = best_scores_[getEntryIdx(best_shift_node, s - 1, m)];
          } else {
            shift_score = -std::numeric_limits<double>::max();
          }

          // a new score is stored as an int, as the node table did; on equal
          // scores the variable path is kept unless it alone has a
          // truncation N-terminal modification
          size_t entry = getEntryIdx(cell, s, m);
          if (var_score > shift_score ||
              (var_score == shift_score && !(var_trunc && !shift_trunc))) {
            if (var_score == -std::numeric_limits<double>::max()) {
              prev_edge_types_[entry] = GRAPH_ALIGN_TYPE_NULL;
              prev_edge_mod_nums_[entry] = 0;
              prev_nodes_[entry] = -1;
              best_scores_[entry] = -std::numeric_limits<int>::max();
            } else {
              prev_edge_types_[entry] = GRAPH_ALIGN_TYPE_VARIABLE;
              prev_edge_mod_nums_[entry] = edge_mod_num;
              prev_nodes_[entry] = best_var_node;
              best_scores_[entry] = toIntScore(var_score + new_score);
              trunc_paths_[entry] = var_trunc;
            }
          } else {
            prev_edge_types_[entry] = GRAPH_ALIGN_TYPE_UNEXPECTED;
            prev_edge_mod_nums_[entry] = 0;
            prev_nodes_[entry] = best_shift_node;
            best_scores_[entry] = toIntScore(shift_score + new_score);
            trunc_paths_[entry] = shift_trunc;
          }
          updateBestShiftNode(i, j, s, m);
        }
      }
    }
  }
  LOG_DEBUG("dp end");
}

GraphResultNodePtrVec GraphAlign::backtrace(int s, int m, bool& trunc) {
  // find the best score; on equal scores, the path without a truncation
  // N-terminal modification
  int best_score = -1;
  int best_node = -1;

  // if whole sequence only
  if (mng_ptr_->whole_protein_only_) {
    size_t cell = getCellIdx(proteo_ver_num_ - 1, spec_ver_num_ - 1);
    int score = toIntScore(best_scores_[getEntryIdx(cell, s, m)]);
    if (score > best_score) {
      best_score = score;
      best_node = static_cast<int>(cell);
    }
  } else {
    for (int i = 0; i < proteo_ver_num_; i++) {
      size_t cell = getCellIdx(i, spec_ver_num_ - 1);
      int score = toIntScore(best_scores_[getEntryIdx(cell, s, m)]);
      if (score > best_score ||
          (score == best_score && isTruncPath(best_node, s, m) &&
           !isTruncPath(static_cast<int>(cell), s, m))) {
        best_score = score;
        best_node = static_cast<int>(cell);
      }
    }
  }
  trunc = best_score > 0 && isTruncPath(best_node, s, m);
  int shift = s;
  int mod = m;
  LOG_DEBUG("best score " << best_score);
  GraphResultNodePtrVec results;
  if (best_score > 0) {
    int cur_node = best_node;
    while (cur_node >= 0) {
      size_t entry = getEntryIdx(cur_node, shift, mod);
      int type = prev_edge_types_[entry];
      results.push_back(std::make_shared<GraphResultNode>(
          cur_node / spec_ver_num_, cur_node % spec_ver_num_, shift, mod,
          type));
      int prev_edge_mod_num = prev_edge_mod_nums_[entry];
      cur_node = prev_nodes_[entry];
      if (type == GRAPH_ALIGN_TYPE_UNEXPECTED) {
        shift--;
      }
      mod = mod - prev_edge_mod_num;
    }
  }
  LOG_DEBUG("obtained result node ptr vec");
  std::reverse(results.begin(), results.end());
  return results;
}

void GraphAlign::backtrace() {
  LOG_DEBUG("start back trace");
  result_nodes_.clear();
  result_truncs_.clear();
  for (int s = 0; s <= mng_ptr_->n_unknown_shift_; s++) {
    LOG_DEBUG("shift num " << s);
    GraphResultNodePtrVec2D nodes_2d;
    std::vector<char> truncs;
    for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
      bool trunc;
      nodes_2d.push_back(backtrace(s, m, trunc));
      truncs.push_back(trunc);
    }
    result_nodes_.push_back(nodes_2d);
    result_truncs_.push_back(truncs);
  }
  LOG_DEBUG("end back trace");
}

void GraphAlign::process() {
  getConsistentPairs();
  initTable();
  dp();
  backtrace();
}

void GraphAlign::getNodeDiagonals(int s, int m) {
  nodes_2d_.clear();
  if (result_nodes_[s][m].size() == 0) {
    return;
  }
  GraphResultNodePtrVec cur_vec;
  cur_vec.push_back(result_nodes_[s][m][0]);
  GraphResultNodePtr prev_node = result_nodes_[s][m][0];
  for (size_t i = 1; i < result_nodes_[s][m].size(); i++) {
    GraphResultNodePtr cur_node = result_nodes_[s][m][i];
    if (cur_node->getShiftNum() == prev_node->getShiftNum() &&
        cur_node->getModNum() == prev_node->getModNum()) {
      cur_vec.push_back(cur_node);
    } else {
      nodes_2d_.push_back(cur_vec);
      cur_vec.clear();
      cur_vec.push_back(cur_node);
    }
    prev_node = cur_node;
  }
  nodes_2d_.push_back(cur_vec);
}

DiagHeaderPtr getFirstDiagonal(const ProteoGraphPtr& proteo_ptr,
                               const GraphResultNodePtrVec& nodes,
                               const std::vector<double>& prm_masses,
                               bool only_diag) {
  int prot_idx = nodes[0]->getFirstIdx();
  int spec_idx = nodes[0]->getSecondIdx();
  double prot_mass = prm_masses[prot_idx];
  double spec_mass = 0;
  if (spec_idx != 0) {
    LOG_ERROR("Spec index is not zero ");
  }
  double shift = spec_mass - prot_mass;
  bool prot_n_term = false;
  bool pep_n_term = false;
  if (prot_idx == 0 || (proteo_ptr->isNme() && prot_idx == 1)) {
    prot_n_term = true;
  } else {
    pep_n_term = true;
  }
  bool prot_c_term = false;
  bool pep_c_term = false;
  int last_node_idx = nodes.size() - 1;
  int last_prot_idx = nodes[last_node_idx]->getFirstIdx();
  if (only_diag) {
    // c_term;
    if (last_prot_idx == static_cast<int>(prm_masses.size()) - 1) {
      prot_c_term = true;
    } else {
      pep_c_term = true;
    }
  }
  DiagHeaderPtr header_ptr = std::make_shared<DiagHeader>(
      shift, true, false, prot_n_term, prot_c_term, pep_n_term, pep_c_term);
  LOG_DEBUG("first diagonal first " << prot_idx << " last " << last_prot_idx);
  header_ptr->setMatchFirstBpPos(prot_idx);
  header_ptr->setMatchLastBpPos(last_prot_idx);
  return header_ptr;
}

DiagHeaderPtr getLastDiagonal(const GraphResultNodePtrVec& nodes,
                              const std::vector<double>& prm_masses,
                              const PrmPeakPtrVec& prm_peaks) {
  int last_node_idx = nodes.size() - 1;
  int last_prot_idx = nodes[last_node_idx]->getFirstIdx();
  int last_spec_idx = nodes[last_node_idx]->getSecondIdx();

  double prot_mass = prm_masses[last_prot_idx];
  double spec_mass = prm_peaks[last_spec_idx]->getPosition();
  double shift = spec_mass - prot_mass;

  bool prot_c_term = false;
  bool pep_c_term = false;
  // c_term;
  if (last_prot_idx == static_cast<int>(prm_masses.size()) - 1) {
    prot_c_term = true;
  } else {
    pep_c_term = true;
  }
  DiagHeaderPtr header_ptr = std::make_shared<DiagHeader>(
      shift, false, true, false, prot_c_term, false, pep_c_term);
  int first_prot_idx = nodes[0]->getFirstIdx();
  LOG_DEBUG("last digaonal first " << first_prot_idx << " last "
                                   << last_prot_idx);
  header_ptr->setMatchFirstBpPos(first_prot_idx);
  header_ptr->setMatchLastBpPos(last_prot_idx);
  return header_ptr;
}

DiagHeaderPtr getInternalDiagonal(const GraphResultNodePtrVec& nodes,
                                  const std::vector<double>& prm_masses,
                                  const PrmPeakPtrVec& prm_peaks) {
  double shift_sum = 0.0;
  for (size_t i = 0; i < nodes.size(); i++) {
    int prot_idx = nodes[i]->getFirstIdx();
    int spec_idx = nodes[i]->getSecondIdx();
    double prot_mass = prm_masses[prot_idx];
    double spec_mass = prm_peaks[spec_idx]->getPosition();
    double shift = spec_mass - prot_mass;
    shift_sum += shift;
  }
  double average_shift = shift_sum / nodes.size();
  DiagHeaderPtr header_ptr = std::make_shared<DiagHeader>(
      average_shift, true, false, false, false, false, false);
  int first_prot_idx = nodes[0]->getFirstIdx();
  header_ptr->setMatchFirstBpPos(first_prot_idx);
  int last_prot_idx = nodes[nodes.size() - 1]->getFirstIdx();
  header_ptr->setMatchLastBpPos(last_prot_idx);
  LOG_DEBUG("internal diagonal first " << first_prot_idx << " last "
                                       << last_prot_idx);
  return header_ptr;
}

void GraphAlign::geneHeaders() {
  diag_headers_.clear();
  diag_headers_2d_.clear();
  std::vector<double> prm_masses =
      proteo_graph_ptr_->getProteoformPtr()->getBpSpecPtr()->getPrmMasses();
  PrmPeakPtrVec prm_peaks = spec_graph_ptr_->getPrmPeakPtrVec();
  if (nodes_2d_.size() >= 1) {
    // add first header
    bool only_diag = false;
    if (nodes_2d_.size() == 1) {
      only_diag = true;
    }
    diag_headers_.push_back(getFirstDiagonal(proteo_graph_ptr_, nodes_2d_[0],
                                             prm_masses, only_diag));
  }
  if (nodes_2d_.size() >= 3) {
    for (size_t i = 1; i < nodes_2d_.size() - 1; i++) {
      diag_headers_.push_back(
          getInternalDiagonal(nodes_2d_[i], prm_masses, prm_peaks));
    }
  }
  if (nodes_2d_.size() >= 2) {
    diag_headers_.push_back(getLastDiagonal(nodes_2d_[nodes_2d_.size() - 1],
                                            prm_masses, prm_peaks));
  }

  // initialize header ptrs
  for (size_t i = 0; i < diag_headers_.size(); i++) {
    double n_shift = diag_headers_[i]->getProtNTermShift();
    double prec_mono_mass =
        spec_graph_ptr_->getSpectrumSetPtr()->getLabelFreePrecMonoMass();
    double c_shift =
        prec_mono_mass -
        proteo_graph_ptr_->getProteoformPtr()->getResSeqPtr()->getSeqMass() -
        n_shift;
    diag_headers_[i]->initHeader(c_shift,
                                 proteo_graph_ptr_->getProteoformPtr());
    LOG_DEBUG("header " << i << " n shift " << n_shift);
  }

  // generate 2d diagonals: first dimension is shift, second dimension is
  // variable mod
  DiagHeaderPtrVec cur_vec;
  for (size_t i = 0; i < diag_headers_.size(); i++) {
    if (nodes_2d_[i][0]->getPrevEdgeType() == GRAPH_ALIGN_TYPE_UNEXPECTED) {
      diag_headers_2d_.push_back(cur_vec);
      cur_vec.clear();
      cur_vec.push_back(diag_headers_[i]);
    } else {
      cur_vec.push_back(diag_headers_[i]);
    }
  }
  diag_headers_2d_.push_back(cur_vec);

  return;
}

PrsmPtr GraphAlign::geneResult(int s, int m) {
  if (result_nodes_[s][m].size() == 0) {
    return nullptr;
  }
  getNodeDiagonals(s, m);
  LOG_DEBUG("begin gene headers ");
  geneHeaders();
  LOG_DEBUG("end gene headers ");
  int last_node_idx = static_cast<int>(result_nodes_[s][m].size() - 1);
  int first_pos = result_nodes_[s][m][0]->getFirstIdx();
  int last_pos = result_nodes_[s][m][last_node_idx]->getFirstIdx() - 1;
  LOG_DEBUG("last pos " << last_pos);
  ProteoformPtr proteo_ptr = proteo_graph_ptr_->getProteoformPtr();

  ProteoformPtr sub_proteo_ptr = toppic::proteoform_factory::geneSubProteoform(
      proteo_ptr, fasta_seq_ptr_, first_pos, last_pos);

  LOG_DEBUG("get sub proteo first pos " << first_pos << " last pos "
                                        << last_pos);
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  ExtendMsPtrVec ms_three_ptr_vec =
      spec_graph_ptr_->getSpectrumSetPtr()->getMsThreePtrVec();
  double min_mass = sp_para_ptr->getMinMass();
  double ppo = sp_para_ptr->getPrecTolerancePtr()->getPpo();
  LOG_DEBUG("begin refine");
  double refine_prec_mass = diagonal_util::refinePrecursorAndHeaderShift(
      proteo_ptr, ms_three_ptr_vec, diag_headers_, ppo, min_mass,
      mng_ptr_->refine_prec_step_width_);
  LOG_DEBUG("get reine prec mass" << refine_prec_mass);

  DeconvMsPtrVec deconv_ms_ptr_vec =
      spec_graph_ptr_->getSpectrumSetPtr()->getDeconvMsPtrVec();
  ExtendMsPtrVec refine_ms_ptr_vec = extend_ms_factory::geneMsThreePtrVec(
      deconv_ms_ptr_vec, sp_para_ptr, refine_prec_mass);

  DiagHeaderPtrVec2D refined_headers_2d = diagonal_util::refineHeadersBgnEnd(
      proteo_ptr, refine_ms_ptr_vec, diag_headers_2d_, diag_headers_, min_mass);

  if (refined_headers_2d.size() == 0) {
    return nullptr;
  }

  DiagHeaderPtrVec refined_headers;

  AlterTypePtrVec shift_types;

  for (size_t i = 0; i < refined_headers_2d.size(); i++) {
    for (size_t j = 0; j < refined_headers_2d[i].size(); j++) {
      refined_headers.push_back(refined_headers_2d[i][j]);
      if (i == 0 && j == 0) {
        shift_types.push_back(nullptr);
      } else if (j == 0) {
        shift_types.push_back(AlterType::UNEXPECTED);
      } else {
        shift_types.push_back(AlterType::VARIABLE);
      }
      LOG_DEBUG("i " << i << " j " << j << " type "
                     << shift_types[shift_types.size() - 1]);
    }
  }

  MassShiftPtrVec shifts = diag_header_util::getDiagonalMassChanges(
      refined_headers, first_pos, last_pos, shift_types);

  sub_proteo_ptr->addMassShiftPtrVec(shifts);
  // sub_proteo_ptr->setVariablePtmNum(m);

  PrsmPtr prsm_ptr = std::make_shared<Prsm>(
      sub_proteo_ptr, deconv_ms_ptr_vec, refine_prec_mass,
      mng_ptr_->prsm_para_ptr_->getSpParaPtr());
  // the graph alignment is TopMG's search, wherever it runs from
  prsm_ptr->setTool("topmg");
  return prsm_ptr;
}

// The PrSM with the most matched fragments over the variable PTM numbers: on
// equal numbers the one with fewer PTMs, unless it alone has a truncation
// N-terminal modification.
PrsmPtr GraphAlign::geneResult(int s) {
  PrsmPtr best_prsm_ptr(nullptr);
  bool best_trunc = false;
  for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
    PrsmPtr cur_prsm_ptr = geneResult(s, m);
    if (cur_prsm_ptr != nullptr) {
      MassShiftPtrVec shift_vec =
          cur_prsm_ptr->getProteoformPtr()->getMassShiftPtrVec(
              AlterType::UNEXPECTED);
      bool valid = true;
      for (size_t i = 0; i < shift_vec.size(); i++) {
        if (std::abs(shift_vec[i]->getMassShift()) > mng_ptr_->max_ptm_mass_) {
          valid = false;
          break;
        }
      }

      bool cur_trunc = result_truncs_[s][m];
      if (valid && (best_prsm_ptr == nullptr ||
                    best_prsm_ptr->getNormMatchTheoFragNum() <
                        cur_prsm_ptr->getNormMatchTheoFragNum() ||
                    (best_prsm_ptr->getNormMatchTheoFragNum() ==
                         cur_prsm_ptr->getNormMatchTheoFragNum() &&
                     best_trunc && !cur_trunc))) {
        best_prsm_ptr = cur_prsm_ptr;
        best_trunc = cur_trunc;
      }
    }
  }
  return best_prsm_ptr;
}

}  // namespace toppic
