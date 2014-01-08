#include <base/logger.hpp>

#include "base/residue_freq.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

ResidueFreq::ResidueFreq(AcidPtrVec acid_list, PtmPtrVec ptm_list,
                         std::string acid_name, std::string ptm_abbr_name, 
                         double freq)
:Residue (acid_list, ptm_list, acid_name, ptm_abbr_name) {

  freq_ = freq;
}

ResFreqPtrVec getResiduePtrVecInstance(AcidPtrVec &acid_list, 
                                           PtmPtrVec &ptm_list,
                                           std::string &file_name) {
  ResFreqPtrVec residue_list;
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int residue_num = getChildCount(parent, "residue_freq");
    LOG_DEBUG( "residue num " << residue_num);
    for (int i = 0; i < residue_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "residue_freq", i);
      std::string acid_name = getChildValue(element, "acid", 0);
      std::string ptm_abbr_name = getChildValue(element, "ptm", 0);
      double freq = getDoubleChildValue(element, "frequency", 0);
      LOG_DEBUG( "acid vec " << acid_list.size() << " ptm vec " 
                << ptm_list.size() << " acid " << acid_name << " ptm " << ptm_abbr_name);
      residue_list.push_back(ResFreqPtr(
              new ResidueFreq(acid_list, ptm_list, acid_name, ptm_abbr_name, freq)));
    }
  }
  return residue_list;
}

}
