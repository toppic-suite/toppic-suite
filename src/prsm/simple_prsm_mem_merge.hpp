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

#ifndef TOPPIC_PRSM_SIMPLE_PRSM_MEM_MERGE_HPP_
#define TOPPIC_PRSM_SIMPLE_PRSM_MEM_MERGE_HPP_

#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace toppic {

// Collects the candidates of the database block threads of a filtering stage
// in memory and writes the four files of the stage, instead of writing one
// file per block and per proteoform type and reading them all back with
// SimplePrsmStrMerge::mergeBlockResults (117 blocks x 4 types for 480_F3_02,
// 87 s single threaded for the variable-PTM filter). This is what commit
// d6f18b30 did for the multiple-shift filter.
//
// The file is byte for byte the one the merge wrote. The merge read the
// records of a spectrum block file by block file, sorted them with
// SimplePrsmStr::cmpScoreDecSeqInc and wrote the first top_num, so the choice
// among records that tie on score and sequence name (the same protein from
// two precursor mass variants) is the one std::sort happens to make on that
// input order. Each candidate therefore carries the position it had in that
// order (its block and the order of the call inside the block), the list is
// restored to it before the final sort, and the same sort on the same records
// makes the same choice.
//
// A spectrum's list is cut back when it grows past twice top_num, but only by
// dropping the records whose score is below the top_num-th score, so no
// record that the merge could have written is lost and the tie behaviour is
// kept. The candidates are held as objects and only the survivors are
// serialized, instead of serializing every candidate of every block.
class SimplePrsmMemMerge {
 public:
  // spec_ids: the spectrum ids of the file in file order
  // (msalign_util::readSpIds). The per-spectrum lists are allocated once, one
  // per spectrum, and never grow while the block threads add to them (a
  // resize would invalidate the references the other threads hold). A
  // spectrum's list is found by its position in the file, not by its id: the
  // ids need not start at 0 or be consecutive, and a list indexed by id up to
  // the number of spectra dropped every candidate of a larger id without a
  // word (issue #598). The lists are written in file order, the order in
  // which the searches read the spectra and the filter file side by side.
  // group_num: the number of candidate groups (the variable PTM filter's
  // rounds), each kept with its own top_nums; the written list of a spectrum
  // is the first group's best, then each later group's best whose protein
  // and range are not listed already (a listed one takes the later
  // candidate's truncation N-terminal modification shifts). One group is
  // the merge above.
  SimplePrsmMemMerge(const std::string& sp_file_name,
                     const std::string& output_file_ext,
                     const std::vector<int>& top_nums,
                     const std::vector<int>& spec_ids, int group_num = 1);

  // type_idx: 0 complete, 1 prefix, 2 suffix, 3 internal
  // block_idx and call_idx give the position the records had in the block
  // files the merge read: block file order, then write order inside a block.
  // group: the candidate group, 0 to group_num - 1
  void add(int type_idx, int block_idx, int call_idx,
           const SimplePrsmPtrVec& prsm_ptrs, int group = 0);

  void write();

 private:
  std::string sp_file_name_;
  std::string output_file_ext_;
  std::vector<int> top_nums_;
  // a candidate with the position it had in the block files
  struct Entry {
    SimplePrsmPtr prsm_ptr_;
    int block_idx_;
    int call_idx_;
    int pos_idx_;
  };

  // the order the candidates of a spectrum are trimmed and written in
  static bool cmpEntry(const Entry& a, const Entry& b);

  // spectrum id -> position of the spectrum in the file; read only once the
  // block threads run
  std::unordered_map<int, int> spec_pos_;
  int group_num_ = 1;
  // results_[group * type number + type][position of the spectrum in the
  // file]
  std::vector<std::vector<std::vector<Entry>>> results_;
  std::vector<std::unique_ptr<std::mutex>> mutexes_;
};

using SimplePrsmMemMergePtr = std::shared_ptr<SimplePrsmMemMerge>;

}  // namespace toppic

#endif
