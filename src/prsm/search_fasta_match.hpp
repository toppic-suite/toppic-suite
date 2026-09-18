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

#ifndef TOPPIC_PRSM_SEARCH_FASTA_MATCH_HPP_
#define TOPPIC_PRSM_SEARCH_FASTA_MATCH_HPP_

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "prsm/prsm.hpp"
#include "seq/fasta_seq.hpp"

namespace toppic {

// Finds every database protein that contains a proteoform's sequence.
//
// The proteins are concatenated into one text and indexed by the positions
// of all their K-mers (K = 4, in a bucket-per-K-mer layout). A query looks up
// the rarest of its K-mers and verifies only those candidate positions, so a
// search costs a few string comparisons instead of a scan of the whole
// database. The result is the same as running std::string::find on every
// protein: the proteins in database order, each with the position of the
// first occurrence of the sequence.
class SearchFastaMatch {
 public:
  explicit SearchFastaMatch(const std::string& db_file_name);

  std::vector<std::pair<FastaSeqPtr, int>> process(const PrsmPtr& prsm_ptr);

  // The proteins containing seq, in database order, with the 0-based
  // position of the first occurrence in each.
  std::vector<std::pair<FastaSeqPtr, int>> findMatches(
      const std::string& seq) const;

 private:
  static constexpr int K = 4;
  // 'A'..'Z' -> 0..25, any other character -> 26
  static constexpr uint32_t ALPHABET_SIZE = 27;

  void buildIndex();

  static uint32_t charCode(char c);

  std::string db_file_name_;
  FastaSeqPtrVec fasta_seq_vec_;
  // all sequences concatenated; seq_start_[i] is where protein i begins and
  // seq_start_[i + 1] where it ends
  std::string text_;
  std::vector<uint32_t> seq_start_;
  // positions_[bucket_start_[c] .. bucket_start_[c + 1]) are the text
  // positions, in increasing order, where the K-mer with code c starts
  std::vector<uint32_t> bucket_start_;
  std::vector<uint32_t> positions_;
};

using SearchFastaMatchPtr = std::shared_ptr<SearchFastaMatch>;

}  // namespace toppic

#endif
