#ifndef PRSM_COMBINE_HPP_
#define PRSM_COMBINE_HPP_

#include <vector>
#include <string>
#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/proteoform_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class PrsmCombine {
 public:
  PrsmCombine(const std::string &db_file_name, 
              const std::string &spec_file_name, 
              const std::vector<std::string> &in_file_exts,
              const std::string &out_file_ext);

  void process();
 private:
  std::string spec_file_name_;
  std::string db_file_name_;
  std::vector<std::string> input_file_exts_;
  std::string output_file_ext_;

};

typedef std::shared_ptr<PrsmCombine> PrsmCombinePtr;
} /* namespace prot */

#endif /* PRSM_COMBINE_HPP_ */
