//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <vector>
#include <map>
#include <string>

#include "base/file_util.hpp"
#include "prsmview/transformer.hpp"
#include "prsmview/anno_view.hpp"

namespace toppic {

void translate(std::map<std::string, std::string> &arguments,
               const std::string &fname_suffix) {
  std::string spectrum_file_name_ = arguments["spectrumFileName"];
  std::string xml_dir = file_util::basename(spectrum_file_name_) + "_" + fname_suffix + "_xml";
  std::string html_dir = file_util::basename(spectrum_file_name_) + "_" + fname_suffix + "_html";
  std::string resource_dir = arguments["resourceDir"];

  file_util::createFolder(html_dir + file_util::getFileSeparator() +"proteoforms");
  file_util::createFolder(html_dir + file_util::getFileSeparator() +"prsms");
  file_util::createFolder(html_dir + file_util::getFileSeparator() +"proteins");
  boost::filesystem::path from_path(resource_dir + file_util::getFileSeparator() + "web");
  boost::filesystem::path to_path(html_dir + file_util::getFileSeparator() + "resources");
  file_util::copyDir(from_path, to_path);

  LOG_DEBUG("trans start!XMLPlatformUtils::Initialize()");
  xercesc::XMLPlatformUtils::Initialize();
  LOG_DEBUG("trans start! XalanTransformer::initialize()");
  xalanc::XalanTransformer::initialize();
  LOG_DEBUG("trans start ! XalanTransformer");
  xalanc::XalanTransformer theXanlanTransformer;

  std::string xml_file_list = xml_dir + file_util::getFileSeparator() + "files.xml";
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

}  // namespace toppic
