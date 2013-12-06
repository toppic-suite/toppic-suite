/*
 * spec_data.cpp
 *
 *  Created on: Dec 5, 2013
 *      Author: xunlikun
 */

#include <log4cxx/logger.h>

#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "spec/spec_data.hpp"

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("SpecData"));

SpecData::SpecData(char const* config_file_name){
	XmlDOMParser* parser = getXmlDOMInstance();
	  if (parser) {
	    LOG4CXX_DEBUG(logger, "config_file_name: " << config_file_name);
	    XmlDOMDocument* doc = new XmlDOMDocument(parser, config_file_name);
	    LOG4CXX_DEBUG(logger, "doc " << doc);
	    if (doc) {
	      xercesc::DOMElement* root = doc->getDocumentElement();
	      LOG4CXX_DEBUG(logger, "root " << root);
	      std::string prm_type_file_name = getChildValue(root, "prm_type_list_file_name", 0);
	      LOG4CXX_DEBUG(logger, "prm type file name: " << prm_type_file_name);
	      prm_peak_type_list_ = getPrmPeakTypePtrVecInstance(prm_type_file_name.c_str());
	      LOG4CXX_DEBUG(logger, "prm type initialized ");

	      std::string support_type_file_name = getChildValue(root, "support_type_list_file_name", 0);
	      LOG4CXX_DEBUG(logger, "support type file name: " << support_type_file_name);
	      support_peak_type_list_ = getSupportPeakTypePtrVecInstance(support_type_file_name.c_str());
	      LOG4CXX_DEBUG(logger, "support type initialized");
	    }
	    delete doc;
	  }
	  // deleting parser is not necessary
}
} /* namespace prot */
