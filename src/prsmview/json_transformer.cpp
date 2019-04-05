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
#include <iostream>
#include <sstream>
#include <fstream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "xml2json/xml2json.hpp"
#include "prsmview/anno_view.hpp"
#include "prsmview/json_transformer.hpp"

namespace toppic {

void jsonConvert(const std::string &xml_file_name, 
                 const std::string &json_file_name) {
  std::ifstream input(xml_file_name);
  std::ostringstream oss;
  oss << input.rdbuf();

  const auto json_str = xml2json( oss.str().data() );
  std::ofstream output(json_file_name);
  output << "prsm_data =" << std::endl;
  output << json_str << std::endl;
  output.close();
}


void jsonTranslate(std::map<std::string, std::string> &arguments,
                   const std::string &fname_suffix) {
  std::string spectrum_file_name_ = arguments["spectrumFileName"];
  std::string xml_dir = file_util::basename(spectrum_file_name_) + "_" + fname_suffix + "_xml";
  std::string html_dir = file_util::basename(spectrum_file_name_) + "_" + fname_suffix + "_html";
  std::string resource_dir = arguments["resourceDir"];

  file_util::createFolder(html_dir + file_util::getFileSeparator() +"proteoforms");
  file_util::createFolder(html_dir + file_util::getFileSeparator() +"prsms");
  file_util::createFolder(html_dir + file_util::getFileSeparator() +"proteins");
  //std::string from_path(resource_dir + file_util::getFileSeparator() + "web");
  //std::string to_path(html_dir + file_util::getFileSeparator() + "resources");
  //file_util::copyDir(from_path, to_path);

  std::string xml_file_list = xml_dir + file_util::getFileSeparator() + "files.xml";
  std::vector<std::vector<std::string>> anno_view = readViewXmlFiles(xml_file_list);

  for (size_t i = 0; i < anno_view.size(); i++) {
    std::cout << "Converting xml files to html files - processing " << i + 1 << " of " << anno_view.size() << " files.\r";
    std::string xml_file_name = anno_view[i][0];
    std::string html_file_name = anno_view[i][2];
    std::string json_file_name = html_file_name.substr(0, html_file_name.length()-4) + "js";

    LOG_DEBUG("xml in " << xml_file_name << " json out " << json_file_name);

    jsonConvert(xml_file_name, json_file_name);
  }

  std::cout << std::endl;
}

}  // namespace toppic

