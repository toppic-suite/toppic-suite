#ifndef PROT_PROTEO_GRAPH_HPP_
#define PROT_PROTEO_GRAPH_HPP_

#include "base/fasta_seq.hpp"
#include "base/proteoform.hpp"
#include "graph/dist_tuple.hpp"
#include "graph/graph.hpp"

namespace prot {

class ProteoGraph {
 public:
  ProteoGraph(FastaSeqPtr seq_ptr, ModPtrVec fix_mod_ptr_vec, MassGraphPtr graph_ptr, 
              bool is_nme, double convert_ratio, int max_mod_num, int max_ptm_sum_mass);

  int getVecIndex(int v1, int v2);
  int getSeqMass(int v1, int v2);

  ProteoformPtr getProteoformPtr() {return db_proteo_ptr_;}
  FastaSeqPtr getFastaSeqPtr() {return db_proteo_ptr_->getFastaSeqPtr();}
  MassGraphPtr getMassGraphPtr() {return graph_ptr_;}
  bool isNme() {return is_nme_;}
  const DistTuplePtrVec2D& getDistTuplePtrVec2D() {return tuple_vec_;}

 private:
  ProteoformPtr db_proteo_ptr_;
  int node_num_;
  int pair_num_;
  std::vector<int> seq_masses_;
  MassGraphPtr graph_ptr_;
  bool is_nme_;
  DistTuplePtrVec2D tuple_vec_;

  void compSeqMasses(double convert_ratio);
  void compDistances(double convert_ratio, int max_mod_num, int max_ptm_sum_mass);
};

typedef std::shared_ptr<ProteoGraph> ProteoGraphPtr;
typedef std::vector<ProteoGraphPtr> ProteoGraphPtrVec;

} /* namespace prot */

#endif /* PROTEO_GRAPH_HPP_ */
