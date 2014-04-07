/*
 * prsm_combine.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: xunlikun
 */

#ifndef PRSM_COMBINE_HPP_
#define PRSM_COMBINE_HPP_

#include <vector>
#include <string>
#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class PrSMCombine {
 public:
  PrSMCombine(std::string db_file, std::string spec_file,
              std::vector<std::string> &in_file_exts, std::string out_file);

  PrSMCombine(std::map<std::string, std::string> arguments,
              std::vector<std::string> &in_file_exts,
              std::string out_file);
  virtual ~PrSMCombine();
  void process();
 private:
  std::string spec_file_;
  std::string db_file_;
  std::vector<std::string> input_file_exts_;
  std::string output_files_;
};
typedef std::shared_ptr<PrSMCombine> PrSMCombinePtr;
} /* namespace prot */

#endif /* PRSM_COMBINE_HPP_ */
