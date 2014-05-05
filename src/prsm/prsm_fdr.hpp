/*
 * prsmfdr.hpp
 *
 *  Created on: Mar 20, 2014
 *      Author: xunlikun
 */

#ifndef PROT_PRSM_FDR_HPP_
#define PROT_PRSM_FDR_HPP_

#include <map>
#include <string>

#include "base/string_util.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class PrsmFdr {
 public:
  PrsmFdr(std::string db_file_name, std::string spec_file_name, 
          std::string input_ext,std::string output_ext);
  PrsmFdr(std::map<std::string,std::string> arguments,std::string input_ext,std::string output_ext);
  void process();
 private:
  std::string db_file_;
  std::string spec_file_;
  std::string input_file_;
  std::string output_file_;

  void compute(PrsmPtrVec &target,PrsmPtrVec decoy);
};

} /* namespace prot */

#endif /* PRSMFDR_HPP_ */
