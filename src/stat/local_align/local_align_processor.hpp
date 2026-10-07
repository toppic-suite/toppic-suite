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

#ifndef TOPPIC_STAT_LOCAL_ALIGN_LOCAL_ALIGN_PROCESSOR_HPP_
#define TOPPIC_STAT_LOCAL_ALIGN_LOCAL_ALIGN_PROCESSOR_HPP_

#include <fstream>
#include <string>

#include "prsm/prsm.hpp"
#include "stat/local/local_mng.hpp"
#include "stat/local_align/local_align_util.hpp"

namespace toppic {

// The alternative to LocalProcessor (toppic --local-align with -B): the
// unexpected shift of a PrSM is explained by the PTMs of the -B file, one
// PTM or an unordered pair, within the precursor tolerance; it is replaced
// by the PTM(s), as variable PTMs with the MIScore of LocalProcessor, only
// when exactly one explanation exists, and without LocalProcessor's filters
// on the matched fragments. A PrSM with two shifts is split at the middle
// matched cleavage site between them and each piece is explained on its
// own.
class LocalAlignProcessor {
 public:
  explicit LocalAlignProcessor(const LocalMngPtr& mng_ptr);

  void process();

 private:
  void init();

  // the switches of TOPPIC_LOCAL_ALIGN_PARA ("name=value,...")
  void readPara();

  PrsmPtr processOnePrsm(const PrsmPtr& prsm_ptr);

  // the shift of the piece explained and localized: the new shifts are
  // appended to new_shifts; false when the shift stays
  bool processPiece(const PrsmPtr& prsm_ptr, const local_align::Piece& piece,
                    const std::string& piece_name, int shift_idx,
                    MassShiftPtrVec& new_shifts);

  void log(const PrsmPtr& prsm_ptr, int shift_idx, double shift_mass,
           const std::string& piece_name, const std::string& expl_names,
           const std::string& decision, const std::string& sites);

  LocalMngPtr mng_ptr_;

  local_align::Para para_;

  std::ofstream log_file_;
};

using LocalAlignProcessorPtr = std::shared_ptr<LocalAlignProcessor>;

}  // namespace toppic

#endif
