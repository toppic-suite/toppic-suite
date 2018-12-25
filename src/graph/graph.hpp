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


#ifndef PROT_GRAPH_HPP_
#define PROT_GRAPH_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "common/base/ptm_base.hpp"
#include "common/base/residue.hpp"
#include "seq/change.hpp"
#include "common/util/logger.hpp"

namespace toppic {

struct EdgeInfo {
  ResiduePtr res_ptr_;

  int change_type_;

  int int_mass_;

  EdgeInfo():res_ptr_(nullptr), change_type_(-1), int_mass_(0) {}

  EdgeInfo(ResiduePtr res_ptr, int change_type, double convert_ratio):
      res_ptr_(res_ptr),
      change_type_(change_type) {
        int_mass_ = static_cast<int>(std::round(res_ptr->getMass() * convert_ratio));
        LOG_DEBUG("int mass " << int_mass_
                  << " res mass " << res_ptr->getMass()
                  << " convert_ratio " << convert_ratio);
      }

  EdgeInfo(double mass, double convert_ratio):res_ptr_(nullptr), change_type_(-1) {
    int_mass_ = static_cast<int>(std::round(mass * convert_ratio));
  }
};

struct VertexInfo {
  int id_;

  VertexInfo():id_(0) {}

  explicit VertexInfo(int id):id_(id) {}
};

struct GraphInfo {
  std::string name_;

  GraphInfo():name_("") {}
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
        VertexInfo, EdgeInfo, GraphInfo> MassGraph;

typedef std::shared_ptr<MassGraph> MassGraphPtr;

typedef std::vector<MassGraphPtr> MassGraphPtrVec;

typedef boost::graph_traits<MassGraph>::vertex_descriptor Vertex;

typedef boost::graph_traits<MassGraph>::edge_descriptor Edge;

}  // namespace toppic

#endif
