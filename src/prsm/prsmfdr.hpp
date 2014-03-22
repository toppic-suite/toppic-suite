/*
 * prsmfdr.hpp
 *
 *  Created on: Mar 20, 2014
 *      Author: xunlikun
 */

#ifndef PRSMFDR_HPP_
#define PRSMFDR_HPP_

#include <map>
#include <string>

#include "base/string_util.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class PrSMFdr {
 public:
  PrSMFdr(std::map<std::string,std::string> arguments,std::string input_ext,std::string output_ext);
  void process();
 private:
  std::string db_file_;
  std::string spec_file_;
  std::string input_file_;
  std::string output_file_;

  void compute(PrSMPtrVec & target,PrSMPtrVec decoy);
};

} /* namespace prot */

#endif /* PRSMFDR_HPP_ */
