// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_GRAPH_DP_NODE_HPP_
#define PROT_GRAPH_DP_NODE_HPP_

#include <memory>
#include <vector>

namespace prot {

#define GRAPH_ALIGN_TYPE_NULL -1
#define GRAPH_ALIGN_TYPE_VARIABLE 0 
#define GRAPH_ALIGN_TYPE_UNEXPECTED 1

class GraphDpNode;
typedef std::shared_ptr<GraphDpNode>  GraphDpNodePtr;
typedef std::vector<GraphDpNodePtr> GraphDpNodePtrVec;
typedef std::vector<GraphDpNodePtrVec> GraphDpNodePtrVec2D;

typedef std::weak_ptr<GraphDpNode>  GraphDpNodeWeakPtr;
typedef std::vector<GraphDpNodeWeakPtr> GraphDpNodeWeakPtrVec;
typedef std::vector<GraphDpNodeWeakPtrVec> GraphDpNodeWeakPtrVec2D;

class GraphDpNode { 
 public:
  GraphDpNode(int first_idx, int second_idx, double node_score,
              int n_unknown_shifts, int max_known_mods);
  int getFirstIdx() {return first_idx_;}
  int getSecondIdx() {return second_idx_;}

  double getNodeScore() {return node_score_;}
  int getPrevEdgeType(int s, int m){return prev_edge_types_[s][m];}
  int getPrevEdgeModNum(int s, int m){return prev_edge_mod_nums_[s][m];}

  double getBestScore(int s, int m) {return best_scores_[s][m];}

  GraphDpNodePtr getPrevNodePtr(int s, int m){return prev_node_ptrs_[s][m].lock();}

  void updateTable(int s, int m, int path_type, int mod_num,
                   GraphDpNodePtr prev_node_ptr, int score);

  void updateBestShiftNode(int s, int m, double score, GraphDpNodePtr prev_node_ptr);

  double getBestShiftScore(int s, int m) {return best_shift_scores_[s][m];}

  GraphDpNodePtr getBestShiftNodePtr(int s, int m) {return best_shift_node_ptrs_[s][m].lock();}


 private:
  int first_idx_;
  int second_idx_;
  // the score for current node
  double node_score_;

  // prev edge and node for dp
  std::vector<std::vector<int>> prev_edge_types_;
  std::vector<std::vector<int>> prev_edge_mod_nums_;

  GraphDpNodeWeakPtrVec2D prev_node_ptrs_;
  std::vector<std::vector<double>> best_scores_;

  // the vector for finding shift nodes
  GraphDpNodeWeakPtrVec2D best_shift_node_ptrs_;
  std::vector<std::vector<double>> best_shift_scores_;
};


} /* namespace prot */

#endif /* GRAPH_DP_NODE_HPP_ */
