#ifndef PROT_PRSM_COVERAGE_HPP_
#define PROT_PRSM_COVERAGE_HPP_

#include "base/string_util.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {

class PrsmConverage {
 public:
  PrsmConverage(PrsmParaPtr prsm_para_ptr, std::string input_file_ext,
                std::string output_file_ext);
  void process();

 private:
  PrsmParaPtr prsm_para_ptr_;
  std::string input_file_ext_;
  std::string output_file_ext_;

};

} /* namespace prot */

#endif /* PROT_PRSM_COVERAGE_HPP_ */
