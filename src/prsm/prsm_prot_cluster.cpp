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

#include "prsm/prsm_prot_cluster.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "common/util/file_util.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace toppic {

namespace prsm_prot_cluster {

void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext) {
  std::string base_name = file_util::basename(spec_file_name);
  std::string input_file_name = base_name + "." + input_file_ext;
  PrsmStrPtrVec prsm_ptrs = prsm_reader_util::readAllPrsmStrs(input_file_name);

  // Group the PrSMs by proteoform cluster; with the PrSMs in E-value order,
  // the first PrSM of a group is the cluster's best proteoform, and the
  // groups come out ordered by that E-value.
  std::sort(prsm_ptrs.begin(), prsm_ptrs.end(), PrsmStr::cmpEValueIncProtInc);
  std::vector<PrsmStrPtrVec> proteo_clusters;
  std::map<int, size_t> proteo_cluster_idx;
  for (const PrsmStrPtr& prsm_ptr : prsm_ptrs) {
    int cluster_id = prsm_ptr->getProteoClusterId();
    auto it = proteo_cluster_idx.find(cluster_id);
    if (it == proteo_cluster_idx.end()) {
      proteo_cluster_idx[cluster_id] = proteo_clusters.size();
      proteo_clusters.push_back(PrsmStrPtrVec{prsm_ptr});
    } else {
      proteo_clusters[it->second].push_back(prsm_ptr);
    }
  }

  // Proteoform clusters whose best proteoforms share a protein form one
  // protein cluster.
  std::map<std::string, int> prot_cluster_ids;
  for (PrsmStrPtrVec& cluster : proteo_clusters) {
    std::string seq_name = cluster[0]->getSeqName();
    auto it = prot_cluster_ids.find(seq_name);
    int prot_cluster_id;
    if (it == prot_cluster_ids.end()) {
      prot_cluster_id = static_cast<int>(prot_cluster_ids.size());
      prot_cluster_ids[seq_name] = prot_cluster_id;
    } else {
      prot_cluster_id = it->second;
    }
    for (PrsmStrPtr& prsm_ptr : cluster) {
      prsm_ptr->setProtClusterId(prot_cluster_id);
    }
  }

  std::sort(prsm_ptrs.begin(), prsm_ptrs.end(),
            PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
  std::string output_file_name = base_name + "." + output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}

}  // namespace prsm_prot_cluster

}  // namespace toppic
