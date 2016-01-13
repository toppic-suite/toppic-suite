#ifndef PROT_GRAPH_HPP_
#define PROT_GRAPH_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "base/ptm_base.hpp"
#include "base/residue.hpp"
#include "base/change.hpp"
#include "base/logger.hpp"

namespace prot {

//#define CONVERT_RATIO 274.335215
//#define CONVERT_RATIO 1

struct EdgeInfo {
  ResiduePtr res_ptr_;
  int change_type_;
  int int_mass_;

  EdgeInfo () {
    res_ptr_ = nullptr;
    change_type_ = -1;
    int_mass_ = 0;
  }
  
  EdgeInfo (ResiduePtr res_ptr, int change_type, double convert_ratio) {
    res_ptr_ = res_ptr;
    change_type_ = change_type;
    int_mass_ = (int)std::round(res_ptr->getMass() * convert_ratio);
    //LOG_DEBUG("int mass " << int_mass_ << " res mass " << res_ptr->getMass() << " convert_ratio " << convert_ratio);
  }

  EdgeInfo (double mass, double convert_ratio) {
    res_ptr_ = nullptr;
    change_type_ = -1;
    int_mass_ = (int)std::round(mass * convert_ratio);
  }
};

struct VertexInfo {
  int id_;

  VertexInfo() {
    id_ = 0;
  }

  VertexInfo(int id) {
    id_ = id;
  }
};

struct GraphInfo {
  std::string name_;

  GraphInfo () {
    name_ = "";
  }
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
        VertexInfo, EdgeInfo, GraphInfo> MassGraph;

typedef std::shared_ptr<MassGraph> MassGraphPtr;

typedef std::vector<MassGraphPtr> MassGraphPtrVec;

typedef boost::graph_traits<MassGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<MassGraph>::edge_descriptor Edge;

template<class T>
class EdgeInfoWriter {
 public:
  EdgeInfoWriter(T _e): e(_e) {}

  template<class VertexOrEdge>
  void operator()(std::ostream& out, const VertexOrEdge& v) const {
    out << "[label=\"";
    if (e[v].res_ptr_ != nullptr) {
      out << e[v].res_ptr_->getAcidPtr()->getOneLetter();
      PtmPtr ptm_ptr = e[v].res_ptr_->getPtmPtr();
      if (!PtmBase::isEmptyPtmPtr(ptm_ptr)) {
        out << "[" << ptm_ptr->getAbbrName() << "]";
      }
      out << ":";
      out << e[v].res_ptr_->getMass() << "\"]";
    }
    else {
      out << e[v].int_mass_ << "\"]";
    }
  }
  private:
    T e;
};

}

#endif
