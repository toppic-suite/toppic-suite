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

// the multiple unexpected shifts stage (-s 2) on the spectra whose first-pass
// PrSM no decoy PrSM could compete with
int TopPIC_multi_ptm_search(std::map<std::string, std::string>& arguments);

// recounting, clustering, PTM characterization, FDR, cutoffs, post mass
// matching and outputs on the best PrSM of every spectrum
// merged_input: the PrSM set to post-process is <base>.toppic_raw_prsm (the
// merged file of a -c run, which has no multiple-shift stage whatever -s
// says) rather than the multiple-shift output.
int TopPIC_multi_ptm_post(std::map<std::string, std::string>& arguments,
                          bool merged_input = false);

int TopPICProgress_multi_file(std::map<std::string, std::string>& arguments,
                              const std::vector<std::string>& spec_file_lst);

}  // namespace toppic

#endif
