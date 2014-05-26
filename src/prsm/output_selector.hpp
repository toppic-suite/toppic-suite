/*
 * output_selector.hpp
 *
 *  Created on: Feb 19, 2014
 *      Author: xunlikun
 */

#ifndef OUTPUT_SELECTOR_HPP_
#define OUTPUT_SELECTOR_HPP_

#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class OutputSelector {
 public:
  OutputSelector(const std::string &db_file,
                 const std::string &spec_file,
                 const std::string &input_file,
                 const std::string &output_file,
                 const std::string &cutoff_type,
                 double cutoff_value);

  OutputSelector(std::map<std::string,std::string> &arguments,
                 const std::string &input_file,
                 const std::string &output_file);
  void process();
 private:
  std::string db_file_;
  std::string spec_file_;
  std::string input_file_;
  std::string output_file_;
  std::string cutoff_type_;
  double cutoff_value_;
};

typedef std::shared_ptr<OutputSelector> OutputSelectorPtr;

} /* namespace prot */

#endif /* OUTPUT_SELECTOR_HPP_ */
