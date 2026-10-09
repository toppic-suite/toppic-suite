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

#include "filter/tag/tag_index.hpp"

#include <algorithm>
#include <cstdlib>

#include "common/base/residue_util.hpp"
#include "common/util/logger.hpp"
#include "seq/fasta_reader.hpp"

namespace toppic {

namespace {

const char kSeparator = '#';

}  // namespace

TagIndex::TagIndex(const std::string& db_file_name,
                   const ModPtrVec& fix_mod_list) {
  FastaReader reader(db_file_name);
  FastaSeqPtr seq_ptr = reader.getNextSeq();
  while (seq_ptr != nullptr) {
    seqs_.push_back(seq_ptr);
    seq_start_.push_back(static_cast<uint32_t>(text_.size()));
    const StringPairVec& pairs = seq_ptr->getAcidPtmPairVec();
    // the residues as the searches build them (proteoform_factory), so the
    // prefix masses are those of the proteoforms the candidates become
    ResiduePtrVec residues =
        residue_util::convertStrToResiduePtrVec(pairs, fix_mod_list);
    if (residues.size() != pairs.size()) {
      LOG_ERROR("The residues of " << seq_ptr->getName()
                                   << " do not match its letters!");
      exit(EXIT_FAILURE);
    }
    double mass = 0;
    for (size_t i = 0; i < pairs.size(); i++) {
      text_.push_back(indexLetter(pairs[i].first[0]));
      prefix_masses_.push_back(mass);
      mass += residues[i]->getMass();
    }
    text_.push_back(kSeparator);
    prefix_masses_.push_back(mass);
    seq_ptr = reader.getNextSeq();
  }
  reader.close();
  seq_start_.push_back(static_cast<uint32_t>(text_.size()));
  buildIndex();
}

uint32_t TagIndex::charCode(char c) {
  c = indexLetter(c);
  if (c >= 'A' && c <= 'Z') {
    return static_cast<uint32_t>(c - 'A');
  }
  return ALPHABET_SIZE - 1;
}

size_t TagIndex::getResidueNum() const {
  // one separator per protein
  return text_.size() - seqs_.size();
}

void TagIndex::buildIndex() {
  size_t bucket_num = 1;
  for (int i = 0; i < K; i++) {
    bucket_num *= ALPHABET_SIZE;
  }
  if (text_.size() < static_cast<size_t>(K)) {
    bucket_start_.assign(bucket_num + 1, 0);
    return;
  }
  size_t window_num = text_.size() - K + 1;
  // the code of every K-mer window, a rolling code: code(p + 1) =
  // (code(p) mod ALPHABET^(K-1)) x ALPHABET + code of the new letter
  uint32_t high = static_cast<uint32_t>(bucket_num / ALPHABET_SIZE);
  std::vector<uint32_t> counts(bucket_num + 1, 0);
  uint32_t code = 0;
  for (int i = 0; i < K; i++) {
    code = code * ALPHABET_SIZE + charCode(text_[i]);
  }
  for (size_t p = 0; p < window_num; p++) {
    counts[code + 1]++;
    if (p + K < text_.size()) {
      code = (code % high) * ALPHABET_SIZE + charCode(text_[p + K]);
    }
  }
  bucket_start_.resize(bucket_num + 1);
  bucket_start_[0] = 0;
  for (size_t c = 0; c < bucket_num; c++) {
    bucket_start_[c + 1] = bucket_start_[c] + counts[c + 1];
  }
  std::vector<uint32_t> next(bucket_start_.begin(), bucket_start_.end() - 1);
  positions_.resize(window_num);
  code = 0;
  for (int i = 0; i < K; i++) {
    code = code * ALPHABET_SIZE + charCode(text_[i]);
  }
  for (size_t p = 0; p < window_num; p++) {
    positions_[next[code]++] = static_cast<uint32_t>(p);
    if (p + K < text_.size()) {
      code = (code % high) * ALPHABET_SIZE + charCode(text_[p + K]);
    }
  }
}

void TagIndex::findTag(const std::string& tag,
                       std::vector<Match>* matches) const {
  matches->clear();
  int len = static_cast<int>(tag.size());
  if (len < K || positions_.empty()) {
    return;
  }
  // the rarest K-mer of the tag
  int best_offset = -1;
  uint32_t best_code = 0;
  size_t best_size = 0;
  for (int offset = 0; offset + K <= len; offset++) {
    uint32_t code = 0;
    for (int i = 0; i < K; i++) {
      code = code * ALPHABET_SIZE + charCode(tag[offset + i]);
    }
    size_t size = bucket_start_[code + 1] - bucket_start_[code];
    if (best_offset < 0 || size < best_size) {
      best_offset = offset;
      best_code = code;
      best_size = size;
    }
  }
  if (best_size == 0) {
    return;
  }
  for (uint32_t k = bucket_start_[best_code]; k < bucket_start_[best_code + 1];
       k++) {
    int64_t start = static_cast<int64_t>(positions_[k]) - best_offset;
    if (start < 0 || start + len > static_cast<int64_t>(text_.size())) {
      continue;
    }
    // a separator in the span cannot match a letter of the tag
    if (text_.compare(static_cast<size_t>(start), len, tag) != 0) {
      continue;
    }
    // the protein that holds the start: the last seq_start_ at or before it
    size_t prot_idx =
        std::upper_bound(seq_start_.begin(), seq_start_.end(),
                         static_cast<uint32_t>(start)) -
        seq_start_.begin() - 1;
    Match match;
    match.prot_idx = static_cast<int>(prot_idx);
    match.pos = static_cast<int>(start - seq_start_[prot_idx]);
    matches->push_back(match);
  }
}

}  // namespace toppic
