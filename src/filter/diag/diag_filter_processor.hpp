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

#ifndef TOPPIC_FILTER_DIAG_DIAG_FILTER_PROCESSOR_HPP_
#define TOPPIC_FILTER_DIAG_DIAG_FILTER_PROCESSOR_HPP_

#include <mutex>
#include <unordered_map>
#include <vector>

#include "filter/diag/diag_filter.hpp"
#include "filter/mng/diag_filter_mng.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace toppic {

class DiagFilterProcessor {
 public:
  DiagFilterProcessor(const DiagFilterMngPtr& mng_ptr) : mng_ptr_(mng_ptr) {};

  void process();

  // Adds the candidates of one spectrum set found in one database block to
  // the running top list of the spectrum (thread-safe).
  void addResults(int spec_id, const SimplePrsmPtrVec& match_ptrs);

  // The spectrum sets converted once for all blocks (empty for topmg's
  // variable-modification variants, which are generated per block).
  const DiagFilterSpectrumPtrVec& getSpectra() const { return spectra_; }

 private:
  DiagFilterMngPtr mng_ptr_;

  DiagFilterSpectrumPtrVec spectra_;

  void readSpectra();

  // spectrum id -> position of the spectrum in the file; read only once the
  // block threads run
  std::unordered_map<int, int> spec_pos_;
  // the best candidates of each spectrum over the blocks processed so far,
  // as the text records that the combined output file receives, by the
  // position of the spectrum in the file: one list per spectrum, allocated
  // before the block threads start and never resized while they hold
  // references to the lists (the ids need not start at 0, issue #598)
  std::vector<SimplePrsmStrPtrVec> spec_results_;
  std::vector<std::mutex> spec_mutexes_;

  void writeResults(const std::string& output_file_name);
};

using DiagFilterProcessorPtr = std::shared_ptr<DiagFilterProcessor>;

}  // namespace toppic

#endif
