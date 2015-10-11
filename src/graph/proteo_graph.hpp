#ifndef PROT_PROTEO_GRAPH_HPP_
#define PROT_PROTEO_GRAPH_HPP_

#include "base/db_residue_seq.hpp"
#include "graph/dist_tuple.hpp"
#include "graph/graph.hpp"

namespace prot {

class ProteoGraph {
 public:
  ProteoGraph(DbResSeqPtr db_res_seq_ptr, MassGraphPtr graph_ptr, 
              bool is_nme, double convert_ratio, int max_ptm_sum_mass);

  int getVecIndex(int v1, int v2);
  int getSeqMass(int v1, int v2);

  ProteoformPtr getProteoformPtr() {return db_proteo_ptr_;}
  DbResSeqPtr getDbResSeqPtr() {return db_proteo_ptr_->getDbResSeqPtr();}
  MassGraphPtr getMassGraphPtr() {return graph_ptr_;}
  bool isNme() {return is_nme_;}
  const DistTuplePtrVec& getDistTuplePtrVec() {return tuple_vec_;}

 private:
  ProteoformPtr db_proteo_ptr_;
  int node_num_;
  int pair_num_;
  std::vector<int> seq_masses_;
  MassGraphPtr graph_ptr_;
  bool is_nme_;
  DistTuplePtrVec tuple_vec_;

  void compSeqMasses(double convert_ratio);
  void compDistances(double convert_ratio, int max_ptm_sum_mass);
};

typedef std::shared_ptr<ProteoGraph> ProteoGraphPtr;
typedef std::vector<ProteoGraphPtr> ProteoGraphPtrVec;

} /* namespace prot */

#endif /* PROTEO_GRAPH_HPP_ */
