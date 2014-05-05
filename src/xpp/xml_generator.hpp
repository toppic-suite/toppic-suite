/*
 * xml_generator.hpp
 *
 *  Created on: Feb 24, 2014
 *      Author: xunlikun
 */

#ifndef XML_GENERATOR_HPP_
#define XML_GENERATOR_HPP_

#include <map>
#include <xercesc/util/PlatformUtils.hpp>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/xml_writer.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "xpp/view_mng.hpp"
#include "xpp/anno_view.hpp"

namespace prot {

class XmlGenerator {
 public:
  XmlGenerator(PrsmParaPtr prsm_para_ptr, std::string exec_dir, 
               std::string input_file);
  void process();
  void outputPrsms(PrSMPtrVec prsms);
  void outputAllPrsms(PrSMPtrVec prsms);
  void outputSpecies(PrSMPtrVec prsms);
  void outputProteins(PrSMPtrVec prsms);
  void outputAllProteins(PrSMPtrVec prsms);
  void outputFileList();

 private:
  std::string input_file_;
  ViewMngPtr mng_;
  ProteoformPtrVec raw_forms_;
  AnnoViewPtr anno_view_;
};

} /* namespace prot */

#endif /* XML_GENERATOR_HPP_ */
