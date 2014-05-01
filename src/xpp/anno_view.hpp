/*
 * anno_view.hpp
 *
 *  Created on: Apr 1, 2014
 *      Author: xunlikun
 */

#ifndef ANNO_VIEW_HPP_
#define ANNO_VIEW_HPP_

#include <map>

#include "base/extreme_value.hpp"
#include "base/proteoform.hpp"
#include "base/anno_residue.hpp"
#include "spec/deconv_peak.hpp"
#include "spec/extend_peak.hpp"
#include "spec/sp_para.hpp"
#include "prsm/cleavage.hpp"
#include "prsm/prsm.hpp"

namespace prot{
class AnnoView {
 public:
  AnnoView();
  std::vector<std::vector<std::string>> file_list_;
  xercesc::DOMElement* geneFileList(XmlDOMDocument* xml_doc);

};

typedef std::shared_ptr<AnnoView> AnnoViewPtr;
std::vector<std::vector<std::string>> readViewXmlFiles(const std::string &file_name);
xercesc::DOMElement* genePrSMView(XmlDOMDocument* xml_doc,PrSMPtr prsm);
xercesc::DOMElement* genePrSMViewAS7(XmlDOMDocument* xml_doc,PrSMPtr prsm);
xercesc::DOMElement* geneProteinView(XmlDOMDocument* xml_doc,ProteoformPtr proteoform,
                                     ExtendMsPtr refine_ms_three,
                                     double min_mass);
xercesc::DOMElement* proteinToXml(XmlDOMDocument* xml_doc,
                                  PrSMPtrVec prsms,
                                  ProteoformPtr protein,
                                  std::vector<int> species);
xercesc::DOMElement* speciesToXml(XmlDOMDocument* xml_doc,PrSMPtrVec prsms);
PrSMPtrVec selectSpeciesPrsms(PrSMPtrVec prsms,int species_id);
std::vector<int> getSpeciesIds(PrSMPtrVec prsms,int seq_id);
std::vector<int> getSpeciesIds(PrSMPtrVec prsms);
xercesc::DOMElement* allProteinToXml(XmlDOMDocument* xml_doc,
                                  PrSMPtrVec prsms,
                                  ProteoformPtrVec proteins);
}

#endif /* ANNO_VIEW_HPP_ */
