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

// Every stage that produces a set of PrSMs returns the file extension the set
// is in, and every later stage is given the extension of the set so far, so
// that TopPICProgress can see which file the post-processing will work on.

// protein filtering + database searching + E-value computation; leaves
// toppic_pre_raw_prsm and returns its extension, or "" if it failed
std::string TopPIC_identify(std::map<std::string, std::string>& arguments);

// TopMG's identification hosted here (hidden --topmg pipeline or
// topmg-only): the ASF filtering with the variable modification file, the
// proteoform-graph alignment, its post-processing, the generating-function
// E-value (the MCMC one with the hidden --topmg-use-mcmc) and the top PrSM
// selection, which chooses between its PrSMs and those of the set it was
// given (none with topmg-only, when raw_prsm_ext is empty and every spectrum
// is searched); returns toppic_topmg_raw_prsm
std::string TopPIC_topmg_identify(
    std::map<std::string, std::string>& arguments,
    const std::string& raw_prsm_ext);

// the multiple unexpected shifts stage (-s 2) on the spectra whose first-pass
// PrSM no decoy PrSM could compete with; merges its PrSMs with the set it was
// given and returns toppic_multi_ptm_raw_prsm
std::string TopPIC_multi_ptm_search(
    std::map<std::string, std::string>& arguments,
    const std::string& raw_prsm_ext);

// recounting, post mass matching, clustering, PTM characterization and FDR on
// the PrSM set of raw_prsm_ext; leaves toppic_record_prsm
int TopPIC_postprocess(std::map<std::string, std::string>& arguments,
                       const std::string& raw_prsm_ext);

// combined runs (-c): clustering and FDR of the merged record PrSMs
// (<name>_ms2.toppic_raw_prsm); leaves <name>_ms2.toppic_record_prsm
int TopPIC_multi_postprocess(std::map<std::string, std::string>& arguments);

// decoy removal (without -K), cutoffs, tables and SQLite output from the
// toppic_record_prsm file
int TopPIC_report(std::map<std::string, std::string>& arguments);

// the two-pass m/z recalibration of one spectrum file (hidden --mass-calib):
// a zero-shift calibration search, the MS1 and MS/MS error models fitted on
// its confident PrSMs (<base>_mass_calib.json, *model_file_name) and the
// corrected spectra <base>_calib_ms2.msalign, whose name it returns (the
// file's own name, uncorrected, with too few calibrants; "" on an error)
std::string TopPIC_mass_calib(std::map<std::string, std::string>& arguments,
                              std::string* model_file_name);

int TopPICProgress_multi_file(std::map<std::string, std::string>& arguments,
                              const std::vector<std::string>& spec_file_lst);

}  // namespace toppic

#endif
