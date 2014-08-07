#ifndef PRSM_SELECTOR_HPP_
#define PRSM_SELECTOR_HPP_

#include <map>

#include "base/string_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class PrsmSelector {
 public:
  PrsmSelector(const std::string &db_file_name,
               const std::string &spec_file_name,
               const std::string &in_file_ext, 
               const std::string &out_file_ext, int n_top);
  void process();
 private:
  std::string spec_file_name_;
  std::string db_file_name_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  int n_top_;

  PrsmPtrVec getTopPrsms(PrsmPtrVec &prsm_ptrs, int n_top);
};

typedef std::shared_ptr<PrsmSelector> PrsmSelectorPtr;

} /* namespace prot */

#endif /* PRSM_SELECTOR_HPP_ */
