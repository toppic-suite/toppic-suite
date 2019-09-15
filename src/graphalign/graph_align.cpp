//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <limits>

#include "seq/proteoform_factory.hpp"
#include "spec/extend_ms_factory.hpp"
#include "oneptmsearch/diagonal.hpp"
#include "oneptmsearch/diagonal_header.hpp"
#include "graph/graph.hpp"
#include "graphalign/graph_align_processor.hpp"
#include "graphalign/graph_align.hpp"

namespace toppic {

std::vector<int> getMinMaxProtDist(DistVec2D dist_vec) {
  std::vector<int> min_dist, max_dist;
  for (size_t i = 0; i < dist_vec.size(); i++) {
    if (dist_vec[i].size() == 0) continue;
    min_dist.push_back(dist_vec[i][0].dist_);
    max_dist.push_back(dist_vec[i][dist_vec[i].size() - 1].dist_);
  }
  std::sort(min_dist.begin(), min_dist.end());
  std::sort(max_dist.begin(), max_dist.end());
  std::vector<int> res;
  res.push_back(min_dist[0]);
  res.push_back(max_dist[max_dist.size() - 1]);
  return res;
}

GraphAlign::GraphAlign(GraphAlignMngPtr mng_ptr,
                       ProteoGraphPtr proteo_graph_ptr,
                       SpecGraphPtr spec_graph_ptr) {
  LOG_DEBUG("Graph constructor start");
  mng_ptr_ = mng_ptr;
  proteo_graph_ptr_ = proteo_graph_ptr;
  spec_graph_ptr_ = spec_graph_ptr;

  dist_vec_ = proteo_graph_ptr_->getDistVec2D();
  spec_dist_ = spec_graph_ptr_->getDistVec();

  std::vector<int> cutoff = getMinMaxProtDist(dist_vec_);
  cutoff[0] -= mng_ptr_->getIntTolerance();
  cutoff[1] += mng_ptr_->getIntTolerance();

  spec_dist_.erase(
      std::remove_if(spec_dist_.begin(), spec_dist_.end(),
                     [cutoff](Dist d){return d.dist_ < cutoff[0] || d.dist_ > cutoff[1];}),
      spec_dist_.end());

  std::sort(spec_dist_.begin(), spec_dist_.end(), distVecUp);

  for (int i = 0; i < mng_ptr->max_known_mods_ + 1; i++) {
    std::sort(dist_vec_[i].begin(), dist_vec_[i].end(), distVecUp);
  }

  pg_ = proteo_graph_ptr_->getMassGraphPtr();
  sg_ = spec_graph_ptr_->getMassGraphPtr();
  proteo_ver_num_ = num_vertices(*pg_.get());
  spec_ver_num_ = num_vertices(*sg_.get());
  LOG_DEBUG("Graph constructor end");
}

void GraphAlign::getConsistentPairs() {
  LOG_DEBUG("consistent pair start");
  int tole = mng_ptr_->getIntTolerance();
  LOG_DEBUG("Integer error tolerance " << tole);
  std::vector<std::pair<int, int>> empty_list;
  std::vector<std::vector<std::pair<int, int>>> empty_vec(mng_ptr_->max_known_mods_ + 1, empty_list);
  for (int i = 0; i < proteo_ver_num_; i++) {
    std::vector<std::vector<std::vector<std::pair<int, int>>>> empty_vec_2d;
    for (int j = 0; j < spec_ver_num_; j++) {
      empty_vec_2d.push_back(empty_vec);
    }
    cons_pairs_.push_back(empty_vec_2d);
  }

  int min_dist = mng_ptr_->getIntMinConsistentDist();
  for (size_t m = 0; m < dist_vec_.size(); m++) {
    if (dist_vec_[m].size() == 0) continue;

    size_t prot_idx_min = 0, prot_idx_max = 0;
    for (size_t spec_idx = 0; spec_idx < spec_dist_.size(); spec_idx++) {
      Dist distance = spec_dist_[spec_idx];
      int sp_dist = distance.dist_;
      if (sp_dist < min_dist) continue;

      bool flag = true;
      while (prot_idx_min < dist_vec_[m].size() && flag) {
        if (dist_vec_[m][prot_idx_min].dist_ >= sp_dist - tole) {
          flag = false;
        } else {
          prot_idx_min++;
        }
      }

      if (prot_idx_min >= dist_vec_[m].size()) continue;

      prot_idx_max = std::max(prot_idx_min, prot_idx_max);

      flag = true;

      while (prot_idx_max < dist_vec_[m].size() && flag) {
        if (dist_vec_[m][prot_idx_max].dist_ > sp_dist + tole) {
          flag = false;
        } else {
          prot_idx_max++;
        }
      }

      for (size_t t = prot_idx_min; t < prot_idx_max; t++) {
        Dist new_distance = dist_vec_[m][t];
        if (std::abs(sp_dist - new_distance.dist_) <= tole)
          addToConsistentPairs(m, distance.pair_ij_, new_distance.pair_ij_);
      }
    }
    dist_vec_[m].clear();
  }
  dist_vec_.clear();
  LOG_DEBUG("consistent pair end");
}

void GraphAlign::addToConsistentPairs(int m, const std::vector<std::pair<int, int>> & sp_pair_ij,
                                      const std::vector<std::pair<int, int>> & pg_pair_ij) {
  for (size_t k = 0; k < pg_pair_ij.size(); k++) {
    for (size_t sp = 0; sp < sp_pair_ij.size(); sp++) {
      int pr_v1 = pg_pair_ij[k].first;
      int pr_v2 = pg_pair_ij[k].second;
      int sp_v1 = sp_pair_ij[sp].first;
      int sp_v2 = sp_pair_ij[sp].second;
      std::pair<int, int> pre_pair(pr_v1, sp_v1);
      cons_pairs_[pr_v2][sp_v2][m].push_back(pre_pair);
    }
  }
}

void GraphAlign::initTable() {
  LOG_DEBUG("init table start");
  double node_score = 1.0;
  for (int i = 0; i < proteo_ver_num_; i++) {
    GraphDpNodePtrVec node_vec;
    for (int j = 0; j < spec_ver_num_; j++) {
      GraphDpNodePtr node_ptr
          = std::make_shared<GraphDpNode>(i, j, node_score, mng_ptr_->n_unknown_shift_,
                                          mng_ptr_->max_known_mods_);
      node_vec.push_back(node_ptr);
    }
    table_.push_back(node_vec);
  }
  LOG_DEBUG("init table step 1");
  for (int i = 0; i < proteo_ver_num_; i++) {
    table_[i][0]->updateTable(0, 0, GRAPH_ALIGN_TYPE_NULL, 0,  nullptr, node_score);
    table_[i][0]->updateBestShiftNode(0, 0, 0, table_[i][0]);
    // LOG_DEBUG("type " << table_[i][0]->getType(0) << " first index " << table_[i][0]->getFirstIdx() << " second index " << table_[i][0]->getSecondIdx());
  }
  LOG_DEBUG("init table end");
}

GraphDpNodePtr GraphAlign::compBestVariableNode(int i, int j, int s, int m, int &best_edge_mod_num) {
  int best_prev_score = -1;
  best_edge_mod_num = -1;
  GraphDpNodePtr best_prev_node = nullptr;
  for (int p = 0; p <= m; p++) {
    for (size_t q = 0; q < cons_pairs_[i][j][p].size(); q++) {
      std::pair<int, int> pair = cons_pairs_[i][j][p][q];
      int pi = pair.first;
      int pj = pair.second;
      // LOG_DEBUG("pi " << pi << " pj " << pj);
      int score = table_[pi][pj]->getBestScore(s, m-p);
      if (score > best_prev_score) {
        best_prev_score = score;
        best_prev_node = table_[pi][pj];
        best_edge_mod_num = p;
      }
    }
  }
  return best_prev_node;
}

GraphDpNodePtr GraphAlign::compBestShiftNode(int i, int j, int s, int m) {
  if (s == 0) {
    return nullptr;
  }
  int best_prev_score = -1;
  GraphDpNodePtr best_prev_node = nullptr;
  GraphDpNodePtr up_node = table_[i-1][j];
  int score = up_node->getBestShiftScore(s-1, m);
  if (score > best_prev_score) {
    best_prev_score = score;
    best_prev_node = up_node->getBestShiftNodePtr(s-1, m);
  }

  GraphDpNodePtr left_node = table_[i][j-1];
  score = left_node->getBestShiftScore(s-1, m);
  if (score > best_prev_score) {
    best_prev_score = score;
    best_prev_node = left_node->getBestShiftNodePtr(s-1, m);
  }
  return best_prev_node;
}

void GraphAlign::updateBestShiftNode(int i, int j, int s, int m) {
  // update best node
  if (table_[i-1][j]->getBestShiftScore(s, m) > table_[i][j-1]->getBestShiftScore(s, m)) {
    table_[i][j]->updateBestShiftNode(s, m, table_[i-1][j]->getBestShiftScore(s, m), table_[i-1][j]->getBestShiftNodePtr(s, m));
  } else {
    table_[i][j]->updateBestShiftNode(s, m, table_[i][j-1]->getBestShiftScore(s, m), table_[i][j-1]->getBestShiftNodePtr(s, m));
  }

  if (table_[i][j]->getBestScore(s, m) > table_[i][j]->getBestShiftScore(s, m)) {
    table_[i][j]->updateBestShiftNode(s, m, table_[i][j]->getBestScore(s, m), table_[i][j]);
  }
}

void GraphAlign::dp() {
  LOG_DEBUG("dp start");
  for (int i = 1; i < proteo_ver_num_; i++) {
    for (int j = 1; j < spec_ver_num_; j++) {
      // compute for zero shift
      for (int s = 0; s <= mng_ptr_->n_unknown_shift_; s++) {
        for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
          int edge_mod_num;
          GraphDpNodePtr best_var_node = compBestVariableNode(i, j, s, m, edge_mod_num);
          double var_score;
          if (best_var_node == nullptr) {
            var_score = - std::numeric_limits<double>::max();
          } else {
            var_score = best_var_node->getBestScore(s, m-edge_mod_num);
          }

          GraphDpNodePtr best_shift_node = compBestShiftNode(i, j, s, m);
          double shift_score;
          if (best_shift_node != nullptr) {
            shift_score = best_shift_node->getBestScore(s-1, m);
          } else {
            shift_score = - std::numeric_limits<double>::max();
          }
          double new_score = table_[i][j]->getNodeScore();
          // LOG_DEBUG("new score " << new_score);

          if (var_score >= shift_score) {
            if (var_score ==  - std::numeric_limits<double>::max()) {
              table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_NULL,
                                        0, nullptr, -std::numeric_limits<int>::max());
            } else {
              table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_VARIABLE,
                                        edge_mod_num, best_var_node, var_score + new_score);
            }
          } else {
            table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_UNEXPECTED,
                                      0, best_shift_node, shift_score + new_score);
          }
          updateBestShiftNode(i, j, s, m);
        }
      }
    }
  }
  LOG_DEBUG("dp end");
}

GraphResultNodePtrVec GraphAlign::backtrace(int s, int m) {
  // find the best score;
  int best_score = -1;
  GraphDpNodePtr best_node_ptr = nullptr;

  for (int i = 0; i < proteo_ver_num_; i++) {
    int score = table_[i][spec_ver_num_-1]->getBestScore(s, m);
    if (score > best_score) {
      best_score = score;
      best_node_ptr = table_[i][spec_ver_num_-1];
    }
  }
  int shift = s;
  int mod = m;
  LOG_DEBUG("best score " << best_score);
  GraphResultNodePtrVec results;
  if (best_score > 0) {
    GraphDpNodePtr cur_node_ptr = best_node_ptr;
    while (cur_node_ptr != nullptr) {
      LOG_DEBUG("cur node " << cur_node_ptr);
      results.push_back(std::make_shared<GraphResultNode>(cur_node_ptr, shift, mod));
      int type = cur_node_ptr->getPrevEdgeType(shift, mod);
      LOG_DEBUG("type " << type << " shift " << shift << " first index " << cur_node_ptr->getFirstIdx() << " second index " << cur_node_ptr->getSecondIdx());
      int prev_edge_mod_num = cur_node_ptr->getPrevEdgeModNum(shift, mod);
      cur_node_ptr = cur_node_ptr->getPrevNodePtr(shift, mod);
      LOG_DEBUG("get prev node ");
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
  for (int s = 0; s <= mng_ptr_->n_unknown_shift_; s++) {
    LOG_DEBUG("shift num " << s);
    GraphResultNodePtrVec2D nodes_2d;
    for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
      nodes_2d.push_back(backtrace(s, m));
    }
    result_nodes_.push_back(nodes_2d);
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
    if (cur_node->getShiftNum() == prev_node->getShiftNum() && cur_node->getModNum() == prev_node->getModNum()) {
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

DiagonalHeaderPtr getFirstDiagonal(ProteoGraphPtr proteo_ptr,
                                   const GraphResultNodePtrVec & nodes,
                                   const std::vector<double> & prm_masses,
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
  DiagonalHeaderPtr header_ptr = std::make_shared<DiagonalHeader>(shift, true, false,
                                                                  prot_n_term, prot_c_term,
                                                                  pep_n_term, pep_c_term);
  LOG_DEBUG("first diagonal first " << prot_idx << " last " << last_prot_idx);
  header_ptr->setMatchFirstBpPos(prot_idx);
  header_ptr->setMatchLastBpPos(last_prot_idx);
  return header_ptr;
}

DiagonalHeaderPtr getLastDiagonal(const GraphResultNodePtrVec & nodes,
                                  const std::vector<double> & prm_masses,
                                  const PrmPeakPtrVec & prm_peaks) {
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
  DiagonalHeaderPtr header_ptr = std::make_shared<DiagonalHeader>(shift, false, true,
                                                                  false, prot_c_term,
                                                                  false, pep_c_term);
  int first_prot_idx = nodes[0]->getFirstIdx();
  LOG_DEBUG("last digaonal first " << first_prot_idx << " last " << last_prot_idx);
  header_ptr->setMatchFirstBpPos(first_prot_idx);
  header_ptr->setMatchLastBpPos(last_prot_idx);
  return header_ptr;
}

DiagonalHeaderPtr getInternalDiagonal(const GraphResultNodePtrVec & nodes,
                                      const std::vector<double> & prm_masses,
                                      const PrmPeakPtrVec & prm_peaks) {
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
  DiagonalHeaderPtr header_ptr
      = std::make_shared<DiagonalHeader>(average_shift, true, false,
                                         false, false, false, false);
  int first_prot_idx = nodes[0]->getFirstIdx();
  header_ptr->setMatchFirstBpPos(first_prot_idx);
  int last_prot_idx = nodes[nodes.size()-1]->getFirstIdx();
  header_ptr->setMatchLastBpPos(last_prot_idx);
  LOG_DEBUG("internal diagonal first " << first_prot_idx << " last " << last_prot_idx);
  return header_ptr;
}

void GraphAlign::geneHeaders() {
  diag_headers_.clear();
  diag_headers_2d_.clear();
  std::vector<double> prm_masses
      = proteo_graph_ptr_->getProteoformPtr()->getBpSpecPtr()->getPrmMasses();
  PrmPeakPtrVec prm_peaks = spec_graph_ptr_->getPrmPeakPtrVec();
  if (nodes_2d_.size() >= 1) {
    // add first header
    bool only_diag = false;
    if (nodes_2d_.size() == 1) {
      only_diag = true;
    }
    diag_headers_.push_back(getFirstDiagonal(proteo_graph_ptr_, nodes_2d_[0], prm_masses, only_diag));
  }
  if (nodes_2d_.size() >= 3) {
    for (size_t i = 1; i < nodes_2d_.size() - 1; i++) {
      diag_headers_.push_back(getInternalDiagonal(nodes_2d_[i], prm_masses, prm_peaks));
    }
  }
  if (nodes_2d_.size() >= 2) {
    diag_headers_.push_back(getLastDiagonal(nodes_2d_[nodes_2d_.size()-1], prm_masses, prm_peaks));
  }

  // initialize header ptrs
  for (size_t i = 0; i < diag_headers_.size(); i++) {
    double n_shift = diag_headers_[i]->getProtNTermShift();
    double prec_mono_mass = spec_graph_ptr_->getSpectrumSetPtr()->getPrecMonoMass();
    double c_shift = prec_mono_mass - proteo_graph_ptr_->getProteoformPtr()->getResSeqPtr()->getSeqMass() - n_shift;
    diag_headers_[i]->initHeader(c_shift,
                                 proteo_graph_ptr_->getProteoformPtr(),
                                 mng_ptr_->align_prefix_suffix_shift_thresh_);
    LOG_DEBUG("header " << i << " n shift " << n_shift);
  }

  // generate 2d diagonals: first dimension is shift, second dimension is
  // variable mod
  DiagonalHeaderPtrVec cur_vec;
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

  ProteoformPtr sub_proteo_ptr
      = toppic::proteoform_factory::geneSubProteoform(proteo_ptr, first_pos, last_pos);

  LOG_DEBUG("get sub proteo first pos " << first_pos << " last pos " << last_pos);
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  ExtendMsPtrVec ms_three_ptr_vec = spec_graph_ptr_->getSpectrumSetPtr()->getMsThreePtrVec();
  double min_mass = sp_para_ptr->getMinMass();
  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  LOG_DEBUG("begin refine");
  double refine_prec_mass = refinePrecursorAndHeaderShift(proteo_ptr, ms_three_ptr_vec,
                                                          diag_headers_, ppo, min_mass,
                                                          mng_ptr_->refine_prec_step_width_);
  LOG_DEBUG("get reine prec mass" << refine_prec_mass);

  DeconvMsPtrVec deconv_ms_ptr_vec = spec_graph_ptr_->getSpectrumSetPtr()->getDeconvMsPtrVec();
  ExtendMsPtrVec refine_ms_ptr_vec
      = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec,  sp_para_ptr, refine_prec_mass);

  DiagonalHeaderPtrVec2D refined_headers_2d = refineHeadersBgnEnd(
      proteo_ptr, refine_ms_ptr_vec, diag_headers_2d_, diag_headers_, min_mass);

  if (refined_headers_2d.size() == 0) {
    return nullptr;
  }

  DiagonalHeaderPtrVec refined_headers;

  AlterTypePtrVec shift_types;

  for (size_t i = 0; i < refined_headers_2d.size(); i++) {
    for (size_t j = 0; j < refined_headers_2d[i].size(); j++) {
      refined_headers.push_back(refined_headers_2d[i][j]);
      if (i == 0 && j == 0) {
        shift_types.push_back(nullptr);
      } else if (j == 0)  {
        shift_types.push_back(AlterType::UNEXPECTED);
      } else {
        shift_types.push_back(AlterType::VARIABLE);
      }
      LOG_DEBUG("i " << i << " j " << j << " type " << shift_types[shift_types.size()-1]);
    }
  }

  MassShiftPtrVec shifts = getDiagonalMassChanges(refined_headers, first_pos, last_pos, shift_types);

  sub_proteo_ptr->addMassShiftPtrVec(shifts);
  sub_proteo_ptr->setVariablePtmNum(m);

  return std::make_shared<Prsm>(sub_proteo_ptr, deconv_ms_ptr_vec, refine_prec_mass,
                                mng_ptr_->prsm_para_ptr_->getSpParaPtr());
}

PrsmPtr GraphAlign::geneResult(int s) {
  PrsmPtr best_prsm_ptr(nullptr);
  for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
    PrsmPtr cur_prsm_ptr = geneResult(s, m);
    if (cur_prsm_ptr != nullptr) {
      MassShiftPtrVec shift_vec
          = cur_prsm_ptr->getProteoformPtr()->getMassShiftPtrVec(AlterType::UNEXPECTED);
      bool valid = true;
      for (size_t i = 0; i < shift_vec.size(); i++) {
        if (std::abs(shift_vec[i]->getMassShift()) > mng_ptr_->max_ptm_mass_) {
          valid = false;
          break;
        }
      }

      if (valid
          && (best_prsm_ptr == nullptr
              || best_prsm_ptr->getNormMatchFragNum() < cur_prsm_ptr->getNormMatchFragNum())) {
        best_prsm_ptr = cur_prsm_ptr;
      }
    }
  }
  return best_prsm_ptr;
}

}  // namespace toppic
