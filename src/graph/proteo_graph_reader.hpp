#ifndef PROT_PROTEO_GRAPH_READER_HPP_
#define PROT_PROTEO_GRAPH_READER_HPP_

#include "graph/graph.hpp"
#include "graph/proteo_anno.hpp"
#include "graph/proteo_graph.hpp"

namespace prot {

class ProteoGraphReader {
 public:
  ProteoGraphReader(const std::string &db_file_name,
                    const ModPtrVec &fix_mod_ptr_vec, 
                    const ProtModPtrVec &prot_mod_ptr_vec,
                    const ModPtrVec &var_mod_ptr_vec,
                    double convert_ratio, int max_mod_num,
                    int max_ptm_sum_mass, int proteo_graph_gap);

  ProteoGraphPtr getNextProteoGraphPtr();

  ProteoGraphPtr getProteoGraphPtrBySeq(FastaSeqPtr seq_ptr);

 private:
  ModPtrVec  fix_mod_ptr_vec_; 
  double convert_ratio_;
  int max_mod_num_;
  int max_ptm_sum_mass_;
  ProteoAnnoPtr proteo_anno_ptr_;
  FastaReaderPtr reader_ptr_;
  int proteo_graph_gap_;

  MassGraphPtr getMassGraphPtr();
};

} /* namespace prot */

#endif /* PROTEO_GRAPH_READER_HPP_ */
