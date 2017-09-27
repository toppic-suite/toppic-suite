// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <vector>

#include "base/file_util.hpp"
#include "prsmview/transformer.hpp"
#include "prsmview/anno_view.hpp"

namespace prot {

void translate(std::map<std::string, std::string> &arguments,
               const std::string &fname_suffix) {
  std::string spectrum_file_name_ = arguments["spectrumFileName"];
  std::string xml_dir = FileUtil::basename(spectrum_file_name_) + "_" + fname_suffix + "_xml";
  std::string html_dir = FileUtil::basename(spectrum_file_name_) + "_" + fname_suffix + "_html";
  std::string exec_dir = arguments["executiveDir"];

  FileUtil::createFolder(html_dir + FileUtil::getFileSeparator() +"proteoforms");
  FileUtil::createFolder(html_dir + FileUtil::getFileSeparator() +"prsms");
  FileUtil::createFolder(html_dir + FileUtil::getFileSeparator() +"proteins");
  boost::filesystem::path from_path(exec_dir + FileUtil::getFileSeparator() + "toppic_resources" + FileUtil::getFileSeparator() + "web");
  boost::filesystem::path to_path(html_dir + FileUtil::getFileSeparator() + "resources");
  FileUtil::copyDir(from_path, to_path);

  LOG_DEBUG("trans start!XMLPlatformUtils::Initialize()");
  xercesc::XMLPlatformUtils::Initialize();
  LOG_DEBUG("trans start! XalanTransformer::initialize()");
  xalanc::XalanTransformer::initialize();
  LOG_DEBUG("trans start ! XalanTransformer");
  xalanc::XalanTransformer theXanlanTransformer;

  std::string xml_file_list = xml_dir + FileUtil::getFileSeparator() + "files.xml";
  std::vector<std::vector<std::string>> anno_view = readViewXmlFiles(xml_file_list);

  for (size_t i = 0; i < anno_view.size(); i++) {
    std::cout << "Converting xml files to html files - processing " << i + 1 << " of " << anno_view.size() << " files.\r";
    const char* xml_in = anno_view[i][0].c_str();
    const char* xsl_in = anno_view[i][1].c_str();
    const char* xml_out = anno_view[i][2].c_str();

    LOG_DEBUG("xml in " << xml_in << " xsl in " << xsl_in << " xml out " << xml_out);

    theXanlanTransformer.transform(xml_in, xsl_in, xml_out);
  }

  std::cout << std::endl;
  xalanc::XalanTransformer::terminate();
  xercesc::XMLPlatformUtils::Terminate();
  xalanc::XalanTransformer::ICUCleanUp();
}

}  // namespace prot
