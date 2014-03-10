/*
 * transformer.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: xunlikun
 */

#include <xpp/transformer.hpp>

Transformer::Transformer() {
  // TODO Auto-generated constructor stub

}

Transformer::~Transformer() {
  // TODO Auto-generated destructor stub
}


void Transformer::trans(){
  xercesc::XMLPlatformUtils::Initialize();
  xalanc::XalanTransformer::initialize();

  xalanc::XalanTransformer theXanlanTransformer;

  const char * xml_in = "xml/prsms.xml";
  const char * xsl_in = "xsl/prsm.xsl";
  const char * xml_out = "html/foo-out.xml";

  int result = theXanlanTransformer.transform(xml_in,xsl_in,xml_out);

  xalanc::XalanTransformer::terminate();
  xercesc::XMLPlatformUtils::Terminate();
  xalanc::XalanTransformer::ICUCleanUp();

}
