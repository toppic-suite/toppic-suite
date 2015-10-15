#ifndef PROT_GRAPH_ALIGN_HPP_
#define PROT_GRAPH_ALIGN_HPP_

#include "prsm/prsm.hpp"
#include "ptmsearch/diagonal_header.hpp"
#include "graph/graph.hpp"
#include "graph/proteo_graph.hpp"
#include "graph/spec_graph.hpp"
#include "graphalign/graph_dp_node.hpp"
#include "graphalign/graph_result_node.hpp"
#include "graphalign/graph_align_mng.hpp"

namespace prot {

typedef std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>> ConsistentPairs;

class GraphAlign {
 public:
  GraphAlign(GraphAlignMngPtr mng_ptr, ProteoGraphPtr proteo_graph_ptr,
             SpecGraphPtr spec_graph_ptr);
  void process();

  PrsmPtr geneResult(int s, int m);
  PrsmPtr geneResult(int s);

 private:
  GraphAlignMngPtr mng_ptr_;
  ProteoGraphPtr proteo_graph_ptr_;
  MassGraphPtr pg_;
  int proteo_ver_num_;
  SpecGraphPtr spec_graph_ptr_;
  MassGraphPtr sg_;
  int spec_ver_num_;
  int max_mod_num_;
  int n_unknown_shift_;
  DistTuplePtrVec2D tuple_vec_;
  ConsistentPairs cons_pairs_;
  GraphDpNodePtrVec2D table_;
  GraphResultNodePtrVec3D result_nodes_;

  GraphResultNodePtrVec2D nodes_2d_;
  DiagonalHeaderPtrVec diag_headers_; 
  DiagonalHeaderPtrVec2D diag_headers_2d_; 

  void getConsistentPairs();
  void initTable();

  GraphDpNodePtr compBestVariableNode(int i, int j, int s, int m, int &best_edge_mod_num);
  GraphDpNodePtr compBestShiftNode(int i, int j, int s, int m);
  void updateBestShiftNode(int i, int j, int s, int m);
  void dp();

  GraphResultNodePtrVec backtrace(int s, int m);
  void backtrace();

  void getNodeDiagonals(int s, int m);

  void geneHeaders();
  
  /*

  //tempary


  */
};

typedef std::shared_ptr<GraphAlign> GraphAlignPtr;

}

#endif

