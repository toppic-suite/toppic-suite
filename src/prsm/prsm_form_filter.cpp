//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include <algorithm>
#include "common/util/file_util.hpp"
#include "seq/fasta_index_reader.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_form_filter.hpp"

namespace toppic {

namespace prsm_form_filter {

void process(const std::string &db_file_name,
             const std::string &spec_file_name,
             const std::string &input_file_ext,
             const std::string &output_file_ext) {

  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;

  ModPtrVec fix_mod_list;
  FastaIndexReaderPtr fasta_reader_ptr 
      = std::make_shared<FastaIndexReader>(db_file_name);
  PrsmPtrVec prsms = prsm_reader_util::readAllPrsms(input_file_name, 
                                                    fasta_reader_ptr,  
                                                    fix_mod_list);

  std::sort(prsms.begin(), prsms.end(), Prsm::cmpEValueIncProtInc);

  PrsmPtrVec selected_forms;

  for (size_t i = 0; i < prsms.size(); i++) {
    // std::cout << "prsm " << i << std::endl;
    bool found = false;
    for (size_t j = 0; j < selected_forms.size(); j++) {
      if (selected_forms[j]->getProteoformPtr()->getProteoClusterId()
          == prsms[i]->getProteoformPtr()->getProteoClusterId()) {
        found = true;
        break;
      }
    }
    if (!found) {
      selected_forms.push_back(prsms[i]);
    }
  }

  // output
  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  std::sort(selected_forms.begin(), selected_forms.end(),
            Prsm::cmpSpecIncPrecIncEvalueIncProtInc);
  writer.writeVector(selected_forms);
  writer.close();
}

}

} /* namespace toppic */
