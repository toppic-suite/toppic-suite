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
//  const xalanc::XalanParsedSource * parsedXml=0;
//  const xalanc::XalanCompiledStylesheet * compliedStyleSheet =0;
//  const xalanc::XSLTInputSource theXSLTInputSource("xsl/proteins.xsl");
//  const xalanc::XSLTInputSource theXMLInputSource("xsl/proteins.xml");
//  xalanc::XSLTInputSource xml_in = "xml/proteins.xml";
//  xalanc::XSLTInputSource xsl_in = "xsl/proteins.xsl";
//  xalanc::XSLTResultTarget xml_out = "html/foo-out.html";
  const char* xml_in = "xml/prsm0.xml";
  const char* xsl_in = "xsl/prsm.xsl";
  const char* xml_out = "html/foo-out.html";

//  theXanlanTransformer.parseSource(theXMLInputSource,parsedXml);
//  theXanlanTransformer.compileStylesheet(theXSLTInputSource,compliedStyleSheet);

  theXanlanTransformer.transform(xml_in,xsl_in,xml_out);

  xalanc::XalanTransformer::terminate();
  xercesc::XMLPlatformUtils::Terminate();
  xalanc::XalanTransformer::ICUCleanUp();

}
