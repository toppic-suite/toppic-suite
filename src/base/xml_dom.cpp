#include <iostream>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <xercesc/framework/XMLFormatter.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMLSSerializer.hpp>
#include <xercesc/dom/DOMLSOutput.hpp>

#include "base/xml_dom.hpp"
#include "base/xml_dom_err_handler.hpp"
 
namespace prot {

XmlDOMParser* XmlDOMParserFactory::dom_parser_ = nullptr;

/* XmlDOMParser */
XmlDOMParser::XmlDOMParser() : parser_(nullptr), err_handler_(nullptr) {
  xercesc::XMLPlatformUtils::Initialize();
  parser_ = new xercesc::XercesDOMParser();
  err_handler_ = (xercesc::ErrorHandler*) new XmlDOMErrorHandler();
  parser_->setErrorHandler(err_handler_);
}

XmlDOMParser::~XmlDOMParser() {
  if (parser_ != nullptr) {
    delete parser_;
    xercesc::XMLPlatformUtils::Terminate();
  }
}

xercesc::DOMDocument* XmlDOMParser::parse(std::string xml_file) {
  parser_->parse(xml_file.c_str());
  return parser_->adoptDocument();
}

/* XmlDOMImplenmation */
XmlDOMImpl* XmlDOMImplFactory::dom_impl_ = nullptr;

XmlDOMImpl::XmlDOMImpl() {
  impl_ = xercesc::DOMImplementationRegistry::getDOMImplementation(X("Core"));
}

XmlDOMImpl::~XmlDOMImpl() {
  if (impl_ != nullptr) {
    delete impl_;
  }
}

xercesc::DOMDocument* XmlDOMImpl::createDoc(std::string root) {
  xercesc::DOMDocument* doc = impl_->createDocument(0, X(root.c_str()), 0);
  return doc;
}

xercesc::DOMLSSerializer* XmlDOMImpl::createSerializer() {
  xercesc::DOMLSSerializer* writer = impl_->createLSSerializer();
  writer->getDomConfig()->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);
  writer->getDomConfig()->setParameter(xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true);
  writer->setNewLine(X("\n"));
  return writer;
}

int writeXmlFile(xercesc::DOMDocument* doc, const char *filename){
  try {
    xercesc::XMLPlatformUtils::Initialize();
  }
  catch(const xercesc::XMLException& e) {
    char* message = xercesc::XMLString::transcode(e.getMessage());
    std::cout << "Error Message: " << message << std::endl;
    xercesc::XMLString::release(&message);
    return 1;
  }
  int errorCode = 0;

  xercesc::DOMImplementation* implementation =  xercesc::DOMImplementationRegistry::getDOMImplementation(X("Core"));

  if (implementation != NULL) {
    try {
      //			xercesc::DOMDocument* doc = impl->createDocument(
      //	                               0,                    // root element namespace URI.
      //	                               X("company"),         // root element name
      //	                               0);                   // document type object (DTD).

      //Return the first registered implementation that has the desired features. In this case, we are after a DOM implementation that has the LS feature... or Load/Save.
      //			xercesc::DOMImplementation *implementation = xercesc::DOMImplementationRegistry::getDOMImplementation(X("LS"));

      // Create a DOMLSSerializer which is used to serialize a DOM tree into an XML document.
      xercesc::DOMLSSerializer *serializer = ((xercesc::DOMImplementationLS*)implementation)->createLSSerializer();

      // Make the output more human readable by inserting line feeds.
      if (serializer->getDomConfig()->canSetParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true)){
        serializer->getDomConfig()->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);
      }

      // The end-of-line sequence of characters to be used in the XML being written out.
      serializer->setNewLine(xercesc::XMLString::transcode("\r\n"));

      // Convert the path into Xerces compatible XMLCh*.
      XMLCh *tempFilePath = xercesc::XMLString::transcode(filename);

      // Specify the target for the XML output.
      xercesc::XMLFormatTarget *formatTarget = new xercesc::LocalFileFormatTarget(tempFilePath);

      // Create a new empty output destination object.
      xercesc::DOMLSOutput *output = ((xercesc::DOMImplementationLS*)implementation)->createLSOutput();

      // Set the stream to our target.
      output->setByteStream(formatTarget);

      // Write the serialized output to the destination.
      serializer->write(doc, output);

      // Cleanup.
      serializer->release();
      xercesc::XMLString::release(&tempFilePath);
      delete formatTarget;
      output->release();

      doc->release();
    }
    catch (const xercesc::OutOfMemoryException&) {
      errorCode = 3;
    }
    catch (const xercesc::DOMException& e) {
      errorCode = 2;
    }
    catch(const xercesc::XMLException& e) {
      char* message = xercesc::XMLString::transcode(e.getMessage());
      std::cout << "Error Message: " << message << std::endl;
      xercesc::XMLString::release(&message);
      return 1;
    }
    catch (...) {
      errorCode = 4;
    }
  }
  xercesc::XMLPlatformUtils::Terminate();
  return errorCode;
}

}
