#include "base/proteoform_factory.hpp"
#include "spec/extend_ms_factory.hpp"
#include "oneptmsearch/diagonal.hpp"
#include "oneptmsearch/diagonal_header.hpp"
#include "graph/graph_util.hpp"
#include "graphalign/graph_align_processor.hpp"
#include "graphalign/graph_align.hpp"

namespace prot {

GraphAlign::GraphAlign(GraphAlignMngPtr mng_ptr, 
                       ProteoGraphPtr proteo_graph_ptr,
                       SpecGraphPtr spec_graph_ptr) {
  LOG_DEBUG("Graph constructor start");
  mng_ptr_ = mng_ptr;
  proteo_graph_ptr_ = proteo_graph_ptr;
  spec_graph_ptr_ = spec_graph_ptr;

  DistTuplePtrVec2D proteo_vec = proteo_graph_ptr_->getDistTuplePtrVec2D();
  tuple_vec_.insert(tuple_vec_.begin(), proteo_vec.begin(), proteo_vec.end());
  /*
  for (size_t i = 0; i < tuple_vec_.size(); i++) {
    LOG_DEBUG("i " << i << " size " << tuple_vec_[i].size());
  }
  for (size_t i = 0; i < tuple_vec_[2].size(); i++) {
    LOG_DEBUG("first " << tuple_vec_[2][i]->getFirstVertex() 
              << " second " << tuple_vec_[2][i]->getSecondVertex()
              << " dist " << tuple_vec_[2][i]->getDist());
  }
  */
  DistTuplePtrVec spec_vec = spec_graph_ptr_->getDistTuplePtrVec();
  for (int i = 0; i < mng_ptr->max_known_mods_ + 1; i++) {
    tuple_vec_[i].insert(tuple_vec_[i].begin(), spec_vec.begin(), spec_vec.end());
    std::sort(tuple_vec_[i].begin(), tuple_vec_[i].end(),distUp);
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
  std::vector<std::pair<int,int>> empty_list;
  std::vector<std::vector<std::pair<int,int>>> empty_vec(mng_ptr_->max_known_mods_ + 1, empty_list);
  for (int i = 0; i < proteo_ver_num_; i++) {
    std::vector<std::vector<std::vector<std::pair<int,int>>>> empty_vec_2d;     
    for (int j = 0; j < spec_ver_num_; j++) {
      empty_vec_2d.push_back(empty_vec);
    }
    cons_pairs_.push_back(empty_vec_2d);
  }

  LOG_DEBUG("tuple vec size " << tuple_vec_.size());
  int min_dist = mng_ptr_->getIntMinConsistentDist();
  for (int m = 0; m < (int)tuple_vec_.size(); m++) {
    for (int i = 0; i < (int)tuple_vec_[m].size(); i++) {
      DistTuplePtr tuple_ptr = tuple_vec_[m][i];
      //LOG_DEBUG("graph_ptr " << tuple_ptr->getGraphPtr());
      if (tuple_ptr->getGraphPtr() == sg_) {
        int sp_v1 = tuple_ptr->getFirstVertex();
        int sp_v2 = tuple_ptr->getSecondVertex();
        int sp_dist = tuple_ptr->getDist();
        // remove zero dist 
        if (sp_dist >= min_dist) {
          //LOG_DEBUG("sp v1 " << sp_v1 << " sp v2 " << sp_v2);
          //search backward
          for (int j = i - 1; j >=0; j--) {
            DistTuplePtr new_tuple_ptr = tuple_vec_[m][j];
            int proteo_dist = new_tuple_ptr->getDist();
            if (sp_dist - proteo_dist > tole) {
              break;
            }
            if (new_tuple_ptr->getGraphPtr() == pg_) {
              int pr_v1 = new_tuple_ptr->getFirstVertex();
              int pr_v2 = new_tuple_ptr->getSecondVertex();
              std::pair<int, int> pre_pair(pr_v1, sp_v1);
              cons_pairs_[pr_v2][sp_v2][m].push_back(pre_pair);
            }
          }
          //search forward
          for (int j = i + 1; j < (int)tuple_vec_[m].size(); j++) {
            DistTuplePtr new_tuple_ptr = tuple_vec_[m][j];
            int proteo_dist = new_tuple_ptr->getDist();
            if (proteo_dist - sp_dist > tole) {
              break;
            }
            if (new_tuple_ptr->getGraphPtr() == pg_) {
              int pr_v1 = new_tuple_ptr->getFirstVertex();
              int pr_v2 = new_tuple_ptr->getSecondVertex();
              //LOG_DEBUG("pr v1 " << pr_v1 << " pr v2 " << pr_v2 << " m " << m << " max mod " << mng_ptr_->max_known_mods_);
              //LOG_DEBUG("m size " << cons_pairs_[pr_v2][sp_v2].size());
              std::pair<int, int> pre_pair(pr_v1, sp_v1);
              //LOG_DEBUG("size " << cons_pairs_[pr_v2][sp_v2][m].size());
              cons_pairs_[pr_v2][sp_v2][m].push_back(pre_pair);
            }
          }
          //LOG_DEBUG("sp v1 sp v2 complete");
        }
      }
    }
  }
  LOG_DEBUG("consistent pair end");
}

void GraphAlign::initTable() {
  LOG_DEBUG("init table start");
  double node_score = 1.0;
  for (int i = 0; i < proteo_ver_num_; i++) {
    GraphDpNodePtrVec node_vec; 
    for (int j = 0; j < spec_ver_num_; j++) {
      GraphDpNodePtr node_ptr(
          new GraphDpNode(i,j,node_score,mng_ptr_->n_unknown_shift_, mng_ptr_->max_known_mods_));
      node_vec.push_back(node_ptr);
    }
    table_.push_back(node_vec);
  }
  LOG_DEBUG("init table step 1");
  for (int i = 0; i < proteo_ver_num_; i++) {
    table_[i][0]->updateTable(0, 0, GRAPH_ALIGN_TYPE_NULL, 0,  nullptr, 0);
    table_[i][0]->updateBestShiftNode(0, 0, 0, table_[i][0]);
    //LOG_DEBUG("type " << table_[i][0]->getType(0) << " first index " << table_[i][0]->getFirstIdx() << " second index " << table_[i][0]->getSecondIdx());
  }
  LOG_DEBUG("init table end");
}

inline GraphDpNodePtr GraphAlign::compBestVariableNode(int i, int j, int s, int m, int &best_edge_mod_num) {
  int best_prev_score = -1;
  best_edge_mod_num = -1;
  GraphDpNodePtr best_prev_node = nullptr;
  for (int p = 0; p <= m; p++) {
    /*
    if (i== 3 && j == 15 && m == 2 && p == 2 ) {
      LOG_DEBUG("consistent pair size "  << cons_pairs_[i][j][p].size());
    }
    */
    for (size_t q = 0; q < cons_pairs_[i][j][p].size(); q++) {
      std::pair<int, int> pair = cons_pairs_[i][j][p][q];
      int pi = pair.first;
      int pj = pair.second;
      //LOG_DEBUG("pi " << pi << " pj " << pj);
      if (table_[pi][pj]->getBestScore(s, m-p) > best_prev_score) {
        best_prev_score = table_[pi][pj]->getBestScore(s, m-p);
        best_prev_node = table_[pi][pj];
        best_edge_mod_num = p;
      }
    }
  }
  return best_prev_node;
}


inline GraphDpNodePtr GraphAlign::compBestShiftNode(int i, int j, int s, int m) {
  if (s == 0) {
    return nullptr;
  }
  int best_prev_score = -1;
  GraphDpNodePtr best_prev_node = nullptr;
  GraphDpNodePtr up_node = table_[i-1][j];
  if (up_node->getBestShiftScore(s-1, m) > best_prev_score) {
    best_prev_score = up_node->getBestShiftScore(s-1, m);
    best_prev_node = up_node->getBestShiftNodePtr(s-1, m);
  }

  GraphDpNodePtr left_node = table_[i][j-1];
  if (left_node->getBestShiftScore(s-1, m) > best_prev_score) {
    best_prev_score = left_node->getBestShiftScore(s-1, m);
    best_prev_node = left_node->getBestShiftNodePtr(s-1, m);
  }
  return best_prev_node;
}

inline void GraphAlign::updateBestShiftNode(int i, int j, int s, int m) {
  // update best node
  if (table_[i-1][j]->getBestShiftScore(s, m) > table_[i][j-1]->getBestShiftScore(s, m)) {
    table_[i][j]->updateBestShiftNode(s, m, table_[i-1][j]->getBestShiftScore(s, m), table_[i-1][j]->getBestShiftNodePtr(s,m));
  }
  else {
    table_[i][j]->updateBestShiftNode(s, m, table_[i][j-1]->getBestShiftScore(s, m), table_[i][j-1]->getBestShiftNodePtr(s,m));
  }
  if (table_[i][j]->getBestScore(s,m) > table_[i][j]->getBestShiftScore(s, m)) {
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

          /*
          if (i == 3) {
            double cur_pos = spec_graph_ptr_->getPrmPeakPtr(j)->getPosition();
            LOG_DEBUG("i " << i << " j " << j << " peak pos " << cur_pos << " s " << s << " m " << m << " var score " << var_score);
          }
          */

          GraphDpNodePtr best_shift_node = compBestShiftNode(i, j, s, m);
          double shift_score;
          if (best_shift_node != nullptr) {
            shift_score = best_shift_node->getBestScore(s-1, m);
          } else {
            shift_score = - std::numeric_limits<double>::max();
          }
          double new_score = table_[i][j]->getNodeScore();
          //LOG_DEBUG("new score " << new_score);

          void updateTable(int s, int m, int path_type, int mod_num,
                           GraphDpNodePtr prev_node_ptr, int score);
          if (var_score >= shift_score) {
            if (var_score ==  - std::numeric_limits<double>::max()) {
              table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_NULL, 
                                        0, nullptr, -std::numeric_limits<double>::max());
            } else {
              table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_VARIABLE,
                                        edge_mod_num, best_var_node, var_score + new_score);
            }
          } 
          else {
            table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_UNEXPECTED,
                                      0, best_shift_node, shift_score + new_score);
          }
          updateBestShiftNode(i,j, s, m);
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
    int score = table_[i][spec_ver_num_-1]->getBestScore(s,m);
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
      results.push_back(GraphResultNodePtr(new GraphResultNode(cur_node_ptr, shift, mod)));
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
  ResSeqPtr res_seq_ptr = proteo_graph_ptr_->getProteoformPtr()->getResSeqPtr();
  for (int i = 0; i < (int)results.size() - 1; i++) {
      int p_cur_idx = results[i]->getFirstIdx();
      int s_cur_idx = results[i]->getSecondIdx();
      int p_prev_idx = results[i+1]->getFirstIdx();
      int s_prev_idx = results[i+1]->getSecondIdx();
      ResSeqPtr sub_seq_ptr = res_seq_ptr->getSubResidueSeq(p_prev_idx, p_cur_idx -1);
      double cur_pos = spec_graph_ptr_->getPrmPeakPtr(s_cur_idx)->getPosition();
      double prev_pos = spec_graph_ptr_->getPrmPeakPtr(s_prev_idx)->getPosition();
      double dist = cur_pos - prev_pos;
      double ptm_mass = dist - sub_seq_ptr->getResMassSum();
      LOG_DEBUG("p cur " << p_cur_idx << " s cur " << s_cur_idx << 
                " p prev " << p_prev_idx << " s prev " << s_prev_idx <<
                " seq " << sub_seq_ptr->toString() << " seq mass " << sub_seq_ptr->getResMassSum() << 
                " cur pos " << cur_pos << " prev pos " << prev_pos << " dist " << dist << 
                " ptm mass " << ptm_mass);
  }
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
    }
    else {
      nodes_2d_.push_back(cur_vec);
      cur_vec.clear();
      cur_vec.push_back(cur_node);
    }
    prev_node = cur_node;
  }
  nodes_2d_.push_back(cur_vec);

  ResSeqPtr res_seq_ptr = proteo_graph_ptr_->getProteoformPtr()->getResSeqPtr();
  for (size_t i = 0; i < nodes_2d_.size(); i++) {
    for (size_t j = 0; j < nodes_2d_[i].size(); j++) {
      int p_idx = nodes_2d_[i][j]->getFirstIdx();
      int s_idx = nodes_2d_[i][j]->getSecondIdx();
      double p_shift = proteo_graph_ptr_->getSeqMass(0, p_idx)/mng_ptr_->convert_ratio_;
      double s_shift = spec_graph_ptr_->getPeakDist(0, s_idx)/mng_ptr_->convert_ratio_;
      double shift = s_shift - p_shift;
      LOG_DEBUG("diag " << i << " index " << j << 
                " p index " << p_idx << " s index " << s_idx <<
                " p shift " << p_shift << " s shift " << s_shift << 
                " N term shift " << shift);
    }
  }
  return; 
}

inline DiagonalHeaderPtr getFirstDiagonal(ProteoGraphPtr proteo_ptr,
                                          GraphResultNodePtrVec nodes, 
                                          std::vector<double> prm_masses, 
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
  }
  else {
    pep_n_term = true;
  }
  bool prot_c_term = false;
  bool pep_c_term = false;
  int last_node_idx = nodes.size() - 1;
  int last_prot_idx = nodes[last_node_idx]->getFirstIdx();
  if (only_diag) {
    // c_term;
    if (last_prot_idx == (int)prm_masses.size() - 1) {
      prot_c_term = true;
    }
    else {
      pep_c_term = true;
    }
  }
  DiagonalHeaderPtr header_ptr(new DiagonalHeader(shift, true, false, 
                                                  prot_n_term, prot_c_term, 
                                                  pep_n_term, pep_c_term));
  LOG_DEBUG("first diagonal first " << prot_idx << " last " << last_prot_idx);
  header_ptr->setMatchFirstBpPos(prot_idx);
  header_ptr->setMatchLastBpPos(last_prot_idx);
  return header_ptr;
}

inline DiagonalHeaderPtr getLastDiagonal(GraphResultNodePtrVec nodes, 
                                         std::vector<double> prm_masses,
                                         PrmPeakPtrVec prm_peaks) {
  int last_node_idx = nodes.size() - 1;
  int last_prot_idx = nodes[last_node_idx]->getFirstIdx();
  int last_spec_idx = nodes[last_node_idx]->getSecondIdx();

  double prot_mass = prm_masses[last_prot_idx];
  double spec_mass = prm_peaks[last_spec_idx]->getPosition();
  double shift = spec_mass - prot_mass;

  bool prot_c_term = false;
  bool pep_c_term = false;
  // c_term;
  if (last_prot_idx == (int)prm_masses.size() - 1) {
    prot_c_term = true;
  }
  else {
    pep_c_term = true;
  }
  DiagonalHeaderPtr header_ptr(new DiagonalHeader(shift, false, true, 
                                                  false, prot_c_term, 
                                                  false, pep_c_term));
  int first_prot_idx = nodes[0]->getFirstIdx();
  LOG_DEBUG("last digaonal first " << first_prot_idx << " last " << last_prot_idx);
  header_ptr->setMatchFirstBpPos(first_prot_idx);
  header_ptr->setMatchLastBpPos(last_prot_idx);
  return header_ptr;
}

inline DiagonalHeaderPtr getInternalDiagonal(GraphResultNodePtrVec nodes, 
                                             std::vector<double> prm_masses,
                                             PrmPeakPtrVec prm_peaks) {
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
  DiagonalHeaderPtr header_ptr(new DiagonalHeader(average_shift, true, false,  
                                                  false, false, false, false));
  int first_prot_idx = nodes[0]->getFirstIdx();
  header_ptr->setMatchFirstBpPos(first_prot_idx);
  int last_prot_idx = nodes[nodes.size()-1]->getFirstIdx();
  header_ptr->setMatchLastBpPos(last_prot_idx);
  LOG_DEBUG("internal diagonal first " << first_prot_idx << " last " << last_prot_idx);
  return header_ptr;
}

void  GraphAlign::geneHeaders() {
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
    for (int i = 1; i < (int)nodes_2d_.size() - 1; i++) {
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
    diag_headers_[i]->initHeader(c_shift, proteo_graph_ptr_->getProteoformPtr(), mng_ptr_->align_prefix_suffix_shift_thresh_);
    LOG_DEBUG("header " << i << " n shift " << n_shift );
  }

  // generate 2d diagonals: first dimension is shift, second dimension is
  // variable mod 
  DiagonalHeaderPtrVec cur_vec;
  for (size_t i = 0; i < diag_headers_.size(); i++) {
    if (nodes_2d_[i][0]->getPrevEdgeType() == GRAPH_ALIGN_TYPE_UNEXPECTED) {
      diag_headers_2d_.push_back(cur_vec);
      cur_vec.clear();
      cur_vec.push_back(diag_headers_[i]);
    }
    else {
      cur_vec.push_back(diag_headers_[i]);
    }
  }
  diag_headers_2d_.push_back(cur_vec);

  return;
}

PrsmPtr GraphAlign::geneResult(int s, int m){
  if(result_nodes_[s][m].size() == 0) {
    return nullptr;
  }
  getNodeDiagonals(s, m);
  LOG_DEBUG("begin gene headers ");
  geneHeaders();
  LOG_DEBUG("end gene headers ");
  int last_node_idx = (int)result_nodes_[s][m].size() - 1;
  int first_pos = result_nodes_[s][m][0]->getFirstIdx();
  int last_pos = result_nodes_[s][m][last_node_idx]->getFirstIdx() - 1;
  LOG_DEBUG("last pos " << last_pos);
  ProteoformPtr proteo_ptr = proteo_graph_ptr_->getProteoformPtr();

  ProteoformPtr sub_proteo_ptr  = ProteoformFactory::geneSubProteoform(proteo_ptr, first_pos, last_pos);

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
      = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec,  sp_para_ptr, refine_prec_mass);
  
  DiagonalHeaderPtrVec2D refined_headers_2d = refineHeadersBgnEnd(
      proteo_ptr, refine_ms_ptr_vec, diag_headers_2d_, diag_headers_, min_mass);

  if(refined_headers_2d.size()==0){
    return nullptr;
  }

  DiagonalHeaderPtrVec refined_headers;
  std::vector<ChangeTypePtr> change_types;
  for (size_t i = 0; i < refined_headers_2d.size(); i++) {
    for (size_t j = 0; j < refined_headers_2d[i].size(); j++) {
      refined_headers.push_back(refined_headers_2d[i][j]);
      if (i == 0 and j == 0) {
        change_types.push_back(nullptr);
      }
      else if (j == 0)  {
        change_types.push_back(ChangeType::UNEXPECTED);
      }
      else {
        change_types.push_back(ChangeType::VARIABLE);
      }
      LOG_DEBUG("i " << i << " j " << j << " type " << change_types[change_types.size()-1]);
    }
  }

  ChangePtrVec changes = getDiagonalMassChanges(refined_headers, first_pos, last_pos, 
                                                change_types);
  /*
  for (size_t i = 0; i < changes.size(); i++) {
    LOG_DEBUG("change " << i << " start " << changes[i]->getLeftBpPos() << " end " << changes[i]->getRightBpPos() << " type " << changes[i]->getChangeType());
  }
  */
  sub_proteo_ptr->addChangePtrVec(changes);

  return PrsmPtr(new Prsm(sub_proteo_ptr, deconv_ms_ptr_vec, refine_prec_mass,
          mng_ptr_->prsm_para_ptr_->getSpParaPtr()));
}

PrsmPtr GraphAlign::geneResult(int s){
  PrsmPtr best_prsm_ptr(nullptr);
  for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
    PrsmPtr cur_prsm_ptr = geneResult(s, m);
    if (cur_prsm_ptr != nullptr) {
      if (best_prsm_ptr == nullptr || best_prsm_ptr->getMatchFragNum() < cur_prsm_ptr->getMatchFragNum()) {
        best_prsm_ptr = cur_prsm_ptr;
      }
    }
  }
  return best_prsm_ptr;
}

}

