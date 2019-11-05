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


#include "seq/proteoform_type.hpp"
#include "prsm/simple_prsm_xml_writer_set.hpp"

namespace toppic {

SimplePrsmXmlWriterSet::SimplePrsmXmlWriterSet(const std::string & output_file_name) { 
  std::string file_name = output_file_name + "_" + ProteoformType::COMPLETE->getName();
  complete_writer_ptr_ = std::make_shared<SimplePrsmXmlWriter>(file_name);
  file_name = output_file_name + "_" + ProteoformType::PREFIX->getName();
  prefix_writer_ptr_ = std::make_shared<SimplePrsmXmlWriter>(file_name);
  file_name = output_file_name + "_" + ProteoformType::SUFFIX->getName();
  suffix_writer_ptr_ = std::make_shared<SimplePrsmXmlWriter>(file_name);
  file_name = output_file_name + "_" + ProteoformType::INTERNAL->getName();
  internal_writer_ptr_ = std::make_shared<SimplePrsmXmlWriter>(file_name);
}

void SimplePrsmXmlWriterSet::close() {
  complete_writer_ptr_->close();
  prefix_writer_ptr_->close();
  suffix_writer_ptr_->close();
  internal_writer_ptr_->close();
}

} /* namespace toppic */
