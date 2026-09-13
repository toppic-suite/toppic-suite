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

#include "prsm/prsm_prot_filter.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <set>
#include <string>

#include "common/base/mod.hpp"
#include "common/util/file_util.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "seq/fasta_index_reader.hpp"

namespace toppic {

namespace prsm_prot_filter {

void process(const std::string& db_file_name, const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext) {
  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;

  ModPtrVec fix_mod_list;
  FastaIndexReaderPtr fasta_reader_ptr =
      std::make_shared<FastaIndexReader>(db_file_name);
  PrsmPtrVec prsms = prsm_reader_util::readAllPrsms(
      input_file_name, fasta_reader_ptr, fix_mod_list);

  // E-value order, so the first PrSM seen for a protein cluster is its best
  // one.
  std::sort(prsms.begin(), prsms.end(), Prsm::cmpEValueIncProtInc);
  PrsmPtrVec selected_prsms;
  std::set<int> seen_prots;
  for (size_t i = 0; i < prsms.size(); i++) {
    int cluster_id = prsms[i]->getProteoformPtr()->getProtClusterId();
    if (seen_prots.insert(cluster_id).second) {
      selected_prsms.push_back(prsms[i]);
    }
  }

  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  std::sort(selected_prsms.begin(), selected_prsms.end(),
            Prsm::cmpSpecIncPrecIncEvalueIncProtInc);
  writer.writeVector(selected_prsms);
  writer.close();
}

}  // namespace prsm_prot_filter

}  // namespace toppic
