#ifndef PROT_PROTEO_GRAPH_READER_HPP_
#define PROT_PROTEO_GRAPH_READER_HPP_

#include "base/db_residue_seq.hpp"
#include "graph/graph.hpp"
#include "graph/proteo_anno.hpp"
#include "graph/proteo_graph.hpp"

namespace prot {

class ProteoGraphReader {
 public:
  ProteoGraphReader(const std::string &db_file_name,
                    const ResiduePtrVec &fix_mod_res_ptr_vec, 
                    const ProtModPtrVec &prot_mod_ptr_vec,
                    const ResiduePtrVec &residue_mod_ptr_vec,
                    double convert_ratio, int max_ptm_sum_mass);

  ProteoGraphPtr getNextProteoGraphPtr();

 private:
  /* db_residue_seq contains fixed modifications */
  ResiduePtrVec fix_mod_res_ptr_vec_; 
  double convert_ratio_;
  int max_ptm_sum_mass_;
  ProteoAnnoPtr proteo_anno_ptr_;
  FastaReaderPtr reader_ptr_;
  int seq_id_ = 0;

  MassGraphPtr getMassGraphPtr();
};

} /* namespace prot */

#endif /* PROTEO_GRAPH_READER_HPP_ */
