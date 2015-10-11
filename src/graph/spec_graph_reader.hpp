#ifndef PROT_SPEC_GRAPH_READER_HPP_
#define PROT_SPEC_GRAPH_READER_HPP_

#include "spec/msalign_reader.hpp"
#include "graph/spec_graph.hpp"

namespace prot {

class SpecGraphReader {
 public:
  SpecGraphReader(const std::string &sp_file_name,
                  int group_sp_num, double convert_ratio,
                  SpParaPtr sp_para_ptr);
  SpecGraphPtr getNextSpecGraphPtr();

 private:
  MsAlignReaderPtr ms_reader_ptr_;
  int group_sp_num_;
  double convert_ratio_;
  SpParaPtr sp_para_ptr_;

  MassGraphPtr getMassGraphPtr(const PrmPeakPtrVec &peak_vec);
};

} /* namespace prot */

#endif /* SPEC_GRAPH_READER_HPP_ */
