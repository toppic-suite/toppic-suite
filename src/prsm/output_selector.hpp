/*
 * output_selector.hpp
 *
 *  Created on: Feb 19, 2014
 *      Author: xunlikun
 */

#ifndef OUTPUT_SELECTOR_HPP_
#define OUTPUT_SELECTOR_HPP_

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class OutputSelector {
 public:
  OutputSelector(std::string db_file,
                 std::string spec_file,
                 std::string input_file,
                 std::string output_file,
                 std::string cutoff_type,
                 double evalue_thresh,
                 double fdr_thresh,
                 double ppo);
  void process();
 private:
  std::string db_file_;
  std::string spec_file_;
  std::string input_file_;
  std::string output_file_;
  std::string cutoff_type_;
  double evalue_thresh_;
  double fdr_thresh_;
  double ppo_;
};

} /* namespace prot */

#endif /* OUTPUT_SELECTOR_HPP_ */
