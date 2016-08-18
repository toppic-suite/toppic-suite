#ifndef PROT_PRSM_PRSM_FORM_FILTER_HPP_
#define PROT_PRSM_PRSM_FORM_FILTER_HPP_

#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace prot {

class PrsmFormFilter {
 public:
  PrsmFormFilter(const std::string &db_file_name,
                 const std::string &spec_file_name,
                 const std::string &input_file_ext,
                 const std::string &output_file_ext,
                 const std::string &output_file_ext_2);
  void process();
 private:
  std::string db_file_name_;
  std::string spec_file_name_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  std::string output_file_ext_2_;
};

typedef std::shared_ptr<PrsmFormFilter> PrsmFormFilterPtr;

} /* namespace prot */

#endif 
