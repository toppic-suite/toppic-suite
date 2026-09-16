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

#ifndef TOPPIC_PRSM_PRSM_PROT_CLUSTER_HPP_
#define TOPPIC_PRSM_PRSM_PROT_CLUSTER_HPP_

#include <string>

namespace toppic {

namespace prsm_prot_cluster {

// Assign the PrSMs of <base>.<input_file_ext> to protein clusters and write
// them to <base>.<output_file_ext>. All PrSMs of a proteoform cluster go to
// the same protein cluster, and two proteoform clusters are merged into one
// protein cluster when their best (lowest E-value) proteoforms are from the
// same protein. Cluster ids are assigned in E-value order of the clusters'
// best PrSMs.
void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext);

}  // namespace prsm_prot_cluster

}  // namespace toppic

#endif
