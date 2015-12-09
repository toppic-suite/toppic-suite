#ifndef PROT_PRSM_PRSM_CUTOFF_SELECTOR_HPP_
#define PROT_PRSM_PRSM_CUTOFF_SELECTOR_HPP_

#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class PrsmCutoffSelector {
 public:
  PrsmCutoffSelector(const std::string &db_file_name,
                     const std::string &spec_file_name,
                     const std::string &input_file_ext,
                     const std::string &output_file_ext,
                     const std::string &cutoff_type,
                     double cutoff_value);
  void process();
 private:
  std::string db_file_name_;
  std::string spec_file_name_;
  std::string cutoff_type_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  double cutoff_value_;
};

typedef std::shared_ptr<PrsmCutoffSelector> PrsmCutoffSelectorPtr;

} /* namespace prot */

#endif /* OUTPUT_SELECTOR_HPP_ */
