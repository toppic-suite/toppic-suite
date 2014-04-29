/*
 * transformer.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: xunlikun
 */

#include <xpp/transformer.hpp>
#include <xpp/anno_view.hpp>
#include "base/command.hpp"

Transformer::Transformer() {
  // TODO Auto-generated constructor stub

}

Transformer::~Transformer() {
  // TODO Auto-generated destructor stub
}


void Transformer::translate(){
	std::cout<<"trans start!XMLPlatformUtils::Initialize()"<<std::endl;
  xercesc::XMLPlatformUtils::Initialize();
  std::cout<<"trans start! XalanTransformer::initialize()"<<std::endl;
  xalanc::XalanTransformer::initialize();
  std::cout<<"trans start ! XalanTransformer"<<std::endl;
  xalanc::XalanTransformer theXanlanTransformer;

  std::string xml_file_list = "xml/files.xml";
  std::vector<std::vector<std::string>> anno_view = prot::readFiles(xml_file_list);
  for(unsigned int i=0;i<anno_view.size();i++){
	  std::cout<<anno_view[i][0]<<std::endl;
    const char* xml_in = anno_view[i][0].c_str();
    const char* xsl_in = anno_view[i][1].c_str();
    const char* xml_out = anno_view[i][2].c_str();

    theXanlanTransformer.transform(xml_in,xsl_in,xml_out);
  }

  xalanc::XalanTransformer::terminate();
  xercesc::XMLPlatformUtils::Terminate();
  xalanc::XalanTransformer::ICUCleanUp();

}
