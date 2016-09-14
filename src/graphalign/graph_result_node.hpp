// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#ifndef PROT_GRAPH_RESULT_NODE_HPP_
#define PROT_GRAPH_RESULT_NODE_HPP_

#include <memory>
#include <vector>

namespace prot {

class GraphResultNode { 
 public:
  GraphResultNode(GraphDpNodePtr node_ptr, int s, int m) {
    first_idx_ = node_ptr->getFirstIdx();
    second_idx_ = node_ptr->getSecondIdx();
    shift_num_ = s;
    mod_num_ = m;
    prev_edge_type_ = node_ptr->getPrevEdgeType(s, m);
  }

  int getFirstIdx() {return first_idx_;}
  int getSecondIdx() {return second_idx_;}
  int getShiftNum(){return shift_num_;}
  int getModNum(){return mod_num_;}
  int getPrevEdgeType() {return prev_edge_type_;}

 private:
  int first_idx_;
  int second_idx_;
  int shift_num_;
  int mod_num_;
  int prev_edge_type_;
};

typedef std::shared_ptr<GraphResultNode> GraphResultNodePtr;
typedef std::vector<GraphResultNodePtr> GraphResultNodePtrVec;
typedef std::vector<GraphResultNodePtrVec> GraphResultNodePtrVec2D;
typedef std::vector<GraphResultNodePtrVec2D> GraphResultNodePtrVec3D;

} /* namespace prot */

#endif /* GRAPH_RESULT_NODE_HPP_ */
