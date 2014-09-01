#ifndef PROT_XML_GENERATOR_HPP_
#define PROT_XML_GENERATOR_HPP_

#include <map>
#include <xercesc/util/PlatformUtils.hpp>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/xml_writer.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "prsmview/prsm_view_mng.hpp"
#include "prsmview/anno_view.hpp"

namespace prot {

class XmlGenerator {
 public:
  XmlGenerator(PrsmParaPtr prsm_para_ptr, std::string exec_dir, 
               std::string input_file);
  void process();
  void outputPrsms(PrsmPtrVec prsms);
  void outputAllPrsms(PrsmPtrVec prsms);
  void outputSpecies(PrsmPtrVec prsms);
  void outputProteins(PrsmPtrVec prsms);
  void outputAllProteins(PrsmPtrVec prsms);
  void outputFileList();

 private:
  std::string input_file_;
  PrsmViewMngPtr mng_;
  ProteoformPtrVec raw_forms_;
  AnnoViewPtr anno_view_;
};

typedef std::shared_ptr<XmlGenerator> XmlGeneratorPtr;

} /* namespace prot */

#endif /* PROT_XML_GENERATOR_HPP_ */
