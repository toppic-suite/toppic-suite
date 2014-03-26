/*
 * xml_generator.hpp
 *
 *  Created on: Feb 24, 2014
 *      Author: xunlikun
 */

#ifndef XML_GENERATOR_HPP_
#define XML_GENERATOR_HPP_

#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/xml_writer.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "xpp/view_mng.hpp"

namespace prot {

class XmlGenerator {
 public:
  XmlGenerator(std::map<std::string,std::string> arguments,std::string input_file);
  void process();
  void processPrSMs(PrSMPtrVec & prsms,ProteoformPtrVec proteoforms);
  void outputPrsms(PrSMPtrVec prsms);
  void outputAllPrsms(PrSMPtrVec prsms);
  void outputProteins(PrSMPtrVec prsms);
  void outputAllProteins(PrSMPtrVec prsms);
 private:
  std::string spec_file_;
  std::string db_file_;
  std::string input_file_;
  std::string output_file_;
  double ppo_;
  ViewMngPtr mng_;
  SpParaPtr sp_para_ptr_;
  std::map<std::string,std::string> arguments_;
};

} /* namespace prot */

#endif /* XML_GENERATOR_HPP_ */
