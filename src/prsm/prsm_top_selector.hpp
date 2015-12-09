#ifndef PROT_PRSM_PRSM_TOP_SELECTOR_HPP_
#define PROT_PRSM_PRSM_TOP_SELECTOR_HPP_

#include <map>

#include "base/string_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace prot {

class PrsmTopSelector {
 public:
  PrsmTopSelector(const std::string &db_file_name,
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

};

typedef std::shared_ptr<PrsmTopSelector> PrsmTopSelectorPtr;

} /* namespace prot */

#endif /* PRSM_SELECTOR_HPP_ */
