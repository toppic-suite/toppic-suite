//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "common/util/str_util.hpp"
#include "prsm/prsm_xml_writer_util.hpp"

namespace toppic {

namespace prsm_xml_writer_util {

PrsmXmlWriterPtrVec geneWriterPtrVec(std::string &file_name, int thread_num) {
  PrsmXmlWriterPtrVec writer_ptr_vec;
  for (int i = 0; i < thread_num; i++) { 
    std::string writer_file_name = file_name + "_" + str_util::toString(i);
    PrsmXmlWriterPtr writer_ptr = std::make_shared<PrsmXmlWriter>(writer_file_name);
    writer_ptr_vec.push_back(writer_ptr);
  }
  return writer_ptr_vec;
}

void closeWriterPtrVec(PrsmXmlWriterPtrVec &writer_ptr_vec) {
  for (size_t i = 0; i < writer_ptr_vec.size(); i++) {
    writer_ptr_vec[i]->close();
  }
}

}

} /* namespace toppic */

