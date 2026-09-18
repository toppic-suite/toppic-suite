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

#include "prsm/search_fasta_match.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "common/util/logger.hpp"
#include "prsm/prsm.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/fasta_seq.hpp"

namespace toppic {

namespace {

// ALPHABET_SIZE^K, the number of K-mer buckets
constexpr uint32_t powAlphabet(int k, uint32_t alphabet_size) {
  return k == 0 ? 1 : alphabet_size * powAlphabet(k - 1, alphabet_size);
}

}  // namespace

SearchFastaMatch::SearchFastaMatch(const std::string& db_file_name)
    : db_file_name_(db_file_name) {
  FastaReader reader(db_file_name_);
  FastaSeqPtr seq_ptr = reader.getNextSeq();
  while (seq_ptr != nullptr) {
    fasta_seq_vec_.push_back(seq_ptr);
    seq_start_.push_back(static_cast<uint32_t>(text_.size()));
    text_ += seq_ptr->getRawSeq();
    seq_ptr = reader.getNextSeq();
  }
  seq_start_.push_back(static_cast<uint32_t>(text_.size()));
  buildIndex();
}

uint32_t SearchFastaMatch::charCode(char c) {
  if (c >= 'A' && c <= 'Z') {
    return static_cast<uint32_t>(c - 'A');
  }
  return ALPHABET_SIZE - 1;
}

void SearchFastaMatch::buildIndex() {
  constexpr uint32_t bucket_num = powAlphabet(K, ALPHABET_SIZE);
  constexpr uint32_t prefix_mod = powAlphabet(K - 1, ALPHABET_SIZE);
  bucket_start_.assign(bucket_num + 1, 0);
  // Visit the K-mers of every protein (never spanning two proteins), first
  // to count the bucket sizes, then to fill the positions in text order.
  for (int pass = 0; pass < 2; pass++) {
    for (size_t i = 0; i + 1 < seq_start_.size(); i++) {
      uint32_t begin = seq_start_[i];
      uint32_t end = seq_start_[i + 1];
      if (end - begin < static_cast<uint32_t>(K)) {
        continue;
      }
      uint32_t code = 0;
      for (uint32_t p = begin; p < end; p++) {
        code = (code % prefix_mod) * ALPHABET_SIZE + charCode(text_[p]);
        if (p + 1 - begin < static_cast<uint32_t>(K)) {
          continue;
        }
        uint32_t start = p + 1 - K;
        if (pass == 0) {
          bucket_start_[code + 1]++;
        } else {
          positions_[bucket_start_[code]++] = start;
        }
      }
    }
    if (pass == 0) {
      for (uint32_t c = 0; c < bucket_num; c++) {
        bucket_start_[c + 1] += bucket_start_[c];
      }
      positions_.resize(bucket_start_[bucket_num]);
    }
  }
  // the second pass advanced every bucket start to its end; shift back
  for (uint32_t c = bucket_num; c > 0; c--) {
    bucket_start_[c] = bucket_start_[c - 1];
  }
  bucket_start_[0] = 0;
  LOG_DEBUG("Indexed " << fasta_seq_vec_.size() << " proteins, "
                       << positions_.size() << " " << K << "-mers");
}

std::vector<std::pair<FastaSeqPtr, int>> SearchFastaMatch::findMatches(
    const std::string& seq) const {
  std::vector<std::pair<FastaSeqPtr, int>> matches;
  size_t m = seq.size();
  if (m == 0) {
    return matches;
  }
  if (m < static_cast<size_t>(K)) {
    // too short for the index: scan each protein
    for (size_t i = 0; i + 1 < seq_start_.size(); i++) {
      std::string_view protein(text_.data() + seq_start_[i],
                               seq_start_[i + 1] - seq_start_[i]);
      size_t found = protein.find(seq);
      if (found != std::string_view::npos) {
        matches.emplace_back(fasta_seq_vec_[i], static_cast<int>(found));
      }
    }
    return matches;
  }
  // the rarest K-mer of the query gives the fewest candidates to verify
  constexpr uint32_t prefix_mod = powAlphabet(K - 1, ALPHABET_SIZE);
  uint32_t code = 0;
  uint32_t best_code = 0;
  size_t best_offset = 0;
  uint32_t best_size = UINT32_MAX;
  for (size_t j = 0; j < m; j++) {
    code = (code % prefix_mod) * ALPHABET_SIZE + charCode(seq[j]);
    if (j + 1 < static_cast<size_t>(K)) {
      continue;
    }
    uint32_t size = bucket_start_[code + 1] - bucket_start_[code];
    if (size < best_size) {
      best_size = size;
      best_code = code;
      best_offset = j + 1 - K;
    }
  }
  if (best_size == 0) {
    return matches;
  }
  // Candidate positions come in increasing text order, so the first hit in
  // each protein is the first occurrence there.
  size_t last_protein = seq_start_.size();
  for (uint32_t b = bucket_start_[best_code]; b < bucket_start_[best_code + 1];
       b++) {
    uint32_t kmer_pos = positions_[b];
    if (kmer_pos < best_offset) {
      continue;
    }
    uint32_t start = kmer_pos - static_cast<uint32_t>(best_offset);
    // the protein containing the K-mer
    size_t i = std::upper_bound(seq_start_.begin(), seq_start_.end(),
                                kmer_pos) -
               seq_start_.begin() - 1;
    if (i == last_protein) {
      continue;
    }
    if (start < seq_start_[i] || start + m > seq_start_[i + 1]) {
      continue;
    }
    if (std::memcmp(text_.data() + start, seq.data(), m) == 0) {
      matches.emplace_back(fasta_seq_vec_[i],
                           static_cast<int>(start - seq_start_[i]));
      last_protein = i;
    }
  }
  return matches;
}

std::vector<std::pair<FastaSeqPtr, int>> SearchFastaMatch::process(
    const PrsmPtr& prsm_ptr) {
  std::string raw_seq =
      prsm_ptr->getProteoformPtr()->getFastaSeqPtr()->getRawSeq();
  int start_pos = prsm_ptr->getProteoformPtr()->getStartPos();
  int end_pos = prsm_ptr->getProteoformPtr()->getEndPos();
  std::string form_seq = raw_seq.substr(start_pos, end_pos - start_pos + 1);
  LOG_DEBUG("proteoform seq " << form_seq);
  return findMatches(form_seq);
}

}  // namespace toppic
