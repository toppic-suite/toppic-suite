//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef XML_BASE_WRITER_HPP_
#define XML_BASE_WRITER_HPP_

#include <iostream>
#include <fstream>
#include <memory>

#include "base/xml_dom_document.hpp"
#include "base/logger.hpp"

namespace prot {

class XmlWriter {
 public:
  XmlWriter(const std::string &file_name,
            const std::string &root);

  ~XmlWriter();

  XmlDOMDocument* getDoc(){return doc_;}

  void write(xercesc::DOMElement* element);

  void write_str(const std::string & str);

  void close();

 private:
  xercesc::DOMLSSerializer * serializer_;

  XmlDOMDocument * doc_;

  std::ofstream file_;

  std::string root_ = "";
};

typedef std::shared_ptr<XmlWriter> XmlWriterPtr;

} /* namespace prot */

#endif /* XML_WRITER_HPP_ */
