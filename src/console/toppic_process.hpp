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

#ifndef TOPPIC_CONSOLE_TOPPIC_PROCESS_HPP_
#define TOPPIC_CONSOLE_TOPPIC_PROCESS_HPP_

#include <map>
#include <string>
#include <vector>

namespace toppic {

int TopPICProgress(std::map<std::string, std::string>& arguments);

// protein filtering + database searching + E-value computation
int TopPIC_identify(std::map<std::string, std::string>& arguments);

// TopMG's identification hosted here (hidden --topmg-identify): the ASF
// filtering with the variable modification file, the proteoform-graph
// alignment, its post-processing, the MCMC E-value and the top PrSM
// selection; its PrSMs are left in toppic_topmg_raw_prsm and are not merged
// into toppic's own results
int TopPIC_topmg_identify(std::map<std::string, std::string>& arguments);

// the multiple unexpected shifts stage (-s 2) on the spectra whose first-pass
// PrSM no decoy PrSM could compete with
int TopPIC_multi_ptm_search(std::map<std::string, std::string>& arguments);

// recounting, post mass matching, clustering, PTM characterization and FDR
// on the best PrSM of every spectrum; leaves toppic_record_prsm
int TopPIC_postprocess(std::map<std::string, std::string>& arguments);

// combined runs (-c): clustering and FDR of the merged record PrSMs
// (<name>_ms2.toppic_raw_prsm); leaves <name>_ms2.toppic_record_prsm
int TopPIC_multi_postprocess(std::map<std::string, std::string>& arguments);

// decoy removal (without -K), cutoffs, tables and SQLite output from the
// toppic_record_prsm file
int TopPIC_report(std::map<std::string, std::string>& arguments);

int TopPICProgress_multi_file(std::map<std::string, std::string>& arguments,
                              const std::vector<std::string>& spec_file_lst);

}  // namespace toppic

#endif
