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

#ifndef TOPPIC_PRSM_PRSM_CLUSTER_PROT_FILTER_HPP_
#define TOPPIC_PRSM_PRSM_CLUSTER_PROT_FILTER_HPP_

#include <string>

namespace toppic {

namespace prsm_cluster_prot_filter {

// Reads the PrSMs of <base_name>.<input_file_ext> (clustered into
// proteoforms, i.e. after prsm_feature_cluster / prsm_simple_cluster), keeps
// in every proteoform cluster only the PrSMs whose protein is the protein of
// the cluster's PrSM with the lowest E-value (ties: the lowest sequence name),
// and writes the rest, in the input order, to <base_name>.<output_file_ext>.
// base_name is the spectrum file name without its extension. Returns the
// number of PrSMs removed.
int process(const std::string& base_name,
            const std::string& input_file_ext,
            const std::string& output_file_ext);

}  // namespace prsm_cluster_prot_filter

}  // namespace toppic

#endif
