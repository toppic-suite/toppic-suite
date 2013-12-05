/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <log4cxx/logger.h>

#include "base/base_data.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("BaseData"));

BaseData::BaseData  (char const* config_file_name) {

  XmlDOMParser* parser = getXmlDOMInstance();
  if (parser) {
    LOG4CXX_DEBUG(logger, "config_file_name: " << config_file_name);
    XmlDOMDocument* doc = new XmlDOMDocument(parser, config_file_name);
    LOG4CXX_DEBUG(logger, "doc " << doc);
    if (doc) {
      xercesc::DOMElement* root = doc->getDocumentElement();
      LOG4CXX_DEBUG(logger, "root " << root);
      std::string acid_file_name = getChildValue(root, "acid_list_file_name", 0);
      LOG4CXX_DEBUG(logger, "acid file name: " << acid_file_name);
      acid_list_ = getAcidPtrVecInstance(acid_file_name.c_str());
      LOG4CXX_DEBUG(logger, "acid initialized ");

      std::string ptm_file_name = getChildValue(root, "ptm_list_file_name", 0);
      LOG4CXX_DEBUG(logger, "ptm file name: " << ptm_file_name);
      ptm_list_ = getPtmPtrVecInstance(ptm_file_name.c_str());
      LOG4CXX_DEBUG(logger, "ptm initialized");

      std::string residue_file_name = getChildValue(root, "residue_list_file_name", 0);
      LOG4CXX_DEBUG(logger, "residue file name: " << residue_file_name);
      residue_list_ = getResiduePtrVecInstance(acid_list_, ptm_list_,
                                               residue_file_name);
      LOG4CXX_DEBUG(logger, "residue initialized");

      std::string trunc_file_name = getChildValue(root, "trunc_list_file_name", 0);
      LOG4CXX_DEBUG(logger, "trunc file name: " << trunc_file_name);
      trunc_list_ = getTruncPtrVecInstance(acid_list_, trunc_file_name.c_str());
      LOG4CXX_DEBUG(logger, "trunc initialized ");

      std::string prot_mod_file_name = getChildValue(root, "prot_mod_list_file_name", 0);
      LOG4CXX_DEBUG(logger, "prot mod file name: " << prot_mod_file_name);
      prot_mod_list_ = getProtModPtrVecInstance(acid_list_, ptm_list_, trunc_list_, prot_mod_file_name.c_str());
      LOG4CXX_DEBUG(logger, "prot mod initialized ");

      std::string ion_type_file_name = getChildValue(root, "ion_type_list_file_name", 0);
      LOG4CXX_DEBUG(logger, "ion type file name: " << ion_type_file_name);
      ion_type_list_ = getIonTypePtrVecInstance(ion_type_file_name.c_str());
      LOG4CXX_DEBUG(logger, "ion type initialized ");

      std::string neutral_loss_file_name = getChildValue(root, "neutral_loss_list_file_name", 0);
      LOG4CXX_DEBUG(logger, "neutral loss file name: " << neutral_loss_file_name);
      neutral_loss_list_ = getNeutralLossPtrVecInstance(neutral_loss_file_name.c_str());
      LOG4CXX_DEBUG(logger, "neutral loss initialized ");

      std::string activation_file_name = getChildValue(root, "activation_list_file_name", 0);
      LOG4CXX_DEBUG(logger, "activation file name: " << activation_file_name);
      activation_list_ = getActivationPtrVecInstance(ion_type_list_, activation_file_name.c_str());
      LOG4CXX_DEBUG(logger, "activation initialized ");

      std::string fix_mod_residue_file_name = getChildValue(root, "fix_mod_residue_file_name", 0);
      LOG4CXX_DEBUG(logger, "fix mod residue file name: " << fix_mod_residue_file_name);
      fix_mod_residue_list_ = getResiduePtrVecInstance(acid_list_, ptm_list_,
                                                       residue_list_, fix_mod_residue_file_name);
      LOG4CXX_DEBUG(logger, "fix mod residue initialized ");

      LOG4CXX_DEBUG(logger, "allow prot mods initialization ");
      xercesc::DOMElement* parent = getChildElement(root, "allow_prot_mod_list", 0);
      int prot_mod_num = getChildCount(parent, "prot_mod");
      for (int i = 0; i < prot_mod_num; i++) {
        std::string mod_name = getChildValue(parent, "prot_mod", i);
        ProtModPtr ptr = getProtModPtrByName(prot_mod_list_, mod_name);
        allow_prot_mod_list_.push_back(ptr);
      }
      LOG4CXX_DEBUG(logger, "allow prot mods initialized ");
    }
    delete doc;
  }
  // deleting parser is not necessary
}

}

