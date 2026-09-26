// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "prsm/prsm_reader.hpp"

#include <cstdlib>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/util/version.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_parser.hpp"
#include "prsm/prsm_util.hpp"

namespace toppic {

PrsmReader::PrsmReader(const std::string& file_name) {
  file_name_ = file_name;
  input_.open(file_name.c_str(), std::ios::in);
}

PrsmReader::~PrsmReader() {
  if (input_.is_open()) {
    input_.close();
  }
}

std::vector<std::string> PrsmReader::readOnePrsmLines() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    str_util::trim(line);
    // LOG_DEBUG("line " << line);
    if (line == "<prsm>") {
      line_list.push_back(line);
    } else if (line == "</prsm>") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    } else if (line == "") {
      continue;
    } else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

PrsmStrPtr PrsmReader::readOnePrsmStr() {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  if (prsm_str_vec.size() == 0) {
    return PrsmStrPtr(nullptr);
  }
  // Result files of TopPIC Suite 1.9 and earlier name the matched numbers
  // match_peak_num / match_fragment_num and have no rt element; PrsmStr
  // would fail on them with an uncaught std::invalid_argument. The tools of
  // this version read only the result files this version writes.
  if (prsm_util::getXmlLine(prsm_str_vec, "<match_exp_frag_num>") == "") {
    LOG_ERROR("The result file "
              << file_name_
              << " was written by an earlier version of TopPIC Suite and "
                 "cannot be read by this version (" << Version::getVersion()
              << "). Run the database search again with this version.");
    exit(EXIT_FAILURE);
  }
  return std::make_shared<PrsmStr>(prsm_str_vec);
}

PrsmPtr PrsmReader::readOnePrsm(const FastaIndexReaderPtr& reader_ptr,
                                const ModPtrVec& fix_mod_list) {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  if (prsm_str_vec.size() == 0) {
    return PrsmPtr(nullptr);
  }
  std::string prsm_str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
  for (size_t i = 0; i < prsm_str_vec.size(); i++) {
    prsm_str += prsm_str_vec[i];
  }
  // LOG_DEBUG("prsm str " << prsm_str);
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  PrsmPtr ptr;
  if (parser) {
    XmlDOMDocument doc(parser->parseStr(prsm_str));
    XmlDOMElement root = doc.getDocumentElement();
    ptr = std::make_shared<Prsm>(root, reader_ptr, fix_mod_list);
  }
  return ptr;
}

void PrsmReader::close() {
  if (input_.is_open()) {
    input_.close();
  }
}

}  // namespace toppic
