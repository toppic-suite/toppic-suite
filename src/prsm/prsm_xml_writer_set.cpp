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
#include "prsm/prsm_xml_writer_set.hpp"

namespace toppic {

PrsmXmlWriterSet::PrsmXmlWriterSet(const std::string & output_file_name, 
                                   int n_mass_shift): 
    n_mass_shift_(n_mass_shift) {
      for (int s = 2; s <= n_mass_shift_; s++) {
        std::string end_str = "_" + str_util::toString(s);
        std::string file_name = output_file_name + "_" + ProteoformType::COMPLETE->getName() + end_str;
        PrsmXmlWriterPtr complete_writer_ptr = std::make_shared<PrsmXmlWriter>(file_name);
        complete_writer_ptrs_.push_back(complete_writer_ptr);
        file_name = output_file_name + "_" + ProteoformType::PREFIX->getName() + end_str;
        PrsmXmlWriterPtr prefix_writer_ptr = std::make_shared<PrsmXmlWriter>(file_name);
        prefix_writer_ptrs_.push_back(prefix_writer_ptr);
        file_name = output_file_name + "_" + ProteoformType::SUFFIX->getName() + end_str;
        PrsmXmlWriterPtr suffix_writer_ptr = std::make_shared<PrsmXmlWriter>(file_name);
        suffix_writer_ptrs_.push_back(suffix_writer_ptr);
        file_name = output_file_name + "_" + ProteoformType::INTERNAL->getName() + end_str;
        PrsmXmlWriterPtr internal_writer_ptr = std::make_shared<PrsmXmlWriter>(file_name);
        internal_writer_ptrs_.push_back(internal_writer_ptr);
      }
    }

void PrsmXmlWriterSet::close() {
  for (int s = 2; s <= n_mass_shift_; s++) {
    complete_writer_ptrs_[s-2]->close();
    prefix_writer_ptrs_[s-2]->close();
    suffix_writer_ptrs_[s-2]->close();
    internal_writer_ptrs_[s-2]->close();
  }
}

} /* namespace toppic */
