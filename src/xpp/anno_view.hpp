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
  std::vector<std::vector<std::string>> file_list_;
  xercesc::DOMElement* geneFileList(XmlDOMDocument* xml_doc);
};

typedef std::shared_ptr<AnnoView> AnnoViewPtr;
std::vector<std::vector<std::string>> readViewXmlFiles(const std::string &file_name);
xercesc::DOMElement* genePrsmView(XmlDOMDocument* xml_doc,PrsmPtr prsm, double min_mass);
xercesc::DOMElement* genePrsmViewAS7(XmlDOMDocument* xml_doc,PrsmPtr prsm);
xercesc::DOMElement* geneProteinView(XmlDOMDocument* xml_doc,ProteoformPtr proteoform,
                                     ExtendMsPtr refine_ms_three, double min_mass);
xercesc::DOMElement* proteinToXml(XmlDOMDocument* xml_doc,
                                  PrsmPtrVec prsms,
                                  ProteoformPtr protein,
                                  std::vector<int> species,
                                  double min_mass);
xercesc::DOMElement* speciesToXml(XmlDOMDocument* xml_doc,PrsmPtrVec prsms, double min_mass);
PrsmPtrVec selectSpeciesPrsms(PrsmPtrVec prsms,int species_id);
std::vector<int> getSpeciesIds(PrsmPtrVec prsms,int seq_id);
std::vector<int> getSpeciesIds(PrsmPtrVec prsms);
xercesc::DOMElement* allProteinToXml(XmlDOMDocument* xml_doc,
                                  PrsmPtrVec prsms,
                                  ProteoformPtrVec proteins,
                                  double min_mass);
}

#endif /* ANNO_VIEW_HPP_ */
