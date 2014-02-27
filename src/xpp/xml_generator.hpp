/*
 * xml_generator.hpp
 *
 *  Created on: Feb 24, 2014
 *      Author: xunlikun
 */

#ifndef XML_GENERATOR_HPP_
#define XML_GENERATOR_HPP_

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/xml_writer.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class XmlGenerator {
 public:
  XmlGenerator(std::string spec_file,std::string db_file,std::string input_file);
  void process();
  void outputPrsms(PrSMPtrVec prsms);
  void outputAllPrsms(PrSMPtrVec prsms);
  void outputProteins(PrSMPtrVec prsms);
  void outputAllProteins(PrSMPtrVec prsms);
 private:
  std::string spec_file_;
  std::string db_file_;
  std::string input_file_;
  std::string output_file_;
};

} /* namespace prot */

#endif /* XML_GENERATOR_HPP_ */
