#ifndef PROT_PRSM_SIMPLE_PRSM_STR_COMBINE_HPP_
#define PROT_PRSM_SIMPLE_PRSM_STR_COMBINE_HPP_

#include <vector>
#include <string>
#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"

namespace prot {

class SimplePrsmStrCombine {
 public:
  SimplePrsmStrCombine(const std::string &spec_file_name, 
                       const std::vector<std::string> &in_file_exts,
                       const std::string &out_file_ext,
                       int top_num);

  SimplePrsmStrCombine(const std::string &spec_file_name, 
                       const std::string &in_file_ext,
                       int in_num,
                       const std::string &out_file_ext, 
                       int top_num);

  void process();
 private:
  std::string spec_file_name_;
  std::vector<std::string> input_file_exts_;
  std::string output_file_ext_;
  int top_num_;
};

typedef std::shared_ptr<SimplePrsmStrCombine> SimplePrsmStrCombinePtr;
} /* namespace prot */

#endif 
