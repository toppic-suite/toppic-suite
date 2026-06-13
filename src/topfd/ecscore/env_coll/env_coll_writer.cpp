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

#include "topfd/ecscore/env_coll/env_coll_writer.hpp"

#include <cstddef>
#include <fstream>
#include <string>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"

namespace toppic {

namespace env_coll_writer {

void writeXmlFeatures(const std::string& output_file_name,
                      const EnvCollPtrVec& env_coll_ptrs) {
  std::ofstream file;
  file.open(output_file_name.c_str());
  LOG_DEBUG("file_name " << output_file_name);
  file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  file << "<env_coll_list>" << std::endl;

  for (size_t i = 0; i < env_coll_ptrs.size(); i++) {
    // pugixml has no detached nodes, so build the env_coll element under a
    // throwaway document root and serialize just that element.
    XmlDOMDocument doc("env_coll_list");
    XmlDOMElement root = doc.getDocumentElement();
    XmlDOMElement element = env_coll_ptrs[i]->toXmlElement(&doc, root);
    std::string str = xml_dom_util::writeToString(element);
    xml_dom_util::writeToStreamByRemovingDoubleLF(file, str);
  }

  file << "</env_coll_list>" << std::endl;
  file.close();
}

}  // namespace env_coll_writer

} /* namespace toppic */
