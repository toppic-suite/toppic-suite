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

#ifndef TOPPIC_FILTER_TAG_TAG_INDEX_HPP_
#define TOPPIC_FILTER_TAG_TAG_INDEX_HPP_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "common/base/mod.hpp"
#include "seq/fasta_seq.hpp"

namespace toppic {

// Every occurrence of a sequence tag in the proteins of a database, and the
// proteins' prefix masses.
//
// The proteins' residue letters (as FastaSeq reads them, B, Z, X and J
// replaced; I written as L, since the two have one mass) are concatenated
// into one text, a separator between proteins, and the positions of every
// K-mer of the text are kept in a bucket per K-mer (prsm/search_fasta_match
// does the same for whole proteoform sequences). A tag of K or more letters
// is looked up by the rarest of its K-mers and verified at each of that
// K-mer's positions, so a lookup costs a few string comparisons. The prefix
// masses (the residues with the fixed modifications, PTMs written in the
// FASTA included) are kept per residue, one double each, so that the
// diagonal shift of a tag occurrence costs one lookup.
class TagIndex {
 public:
  static constexpr int K = 4;

  // an occurrence: the protein (its index in the database file) and the
  // 0-based position of the tag's first residue in it
  struct Match {
    int prot_idx;
    int pos;
  };

  // db_file_name: the search database (targets and decoys), a FASTA file;
  // fix_mod_list: the fixed modifications of the residue masses
  TagIndex(const std::string& db_file_name, const ModPtrVec& fix_mod_list);

  // The occurrences of tag (K or more letters, I written as L), in the
  // order of the text. Letters outside A-Z never match.
  void findTag(const std::string& tag, std::vector<Match>* matches) const;

  int getProtNum() const { return static_cast<int>(seqs_.size()); }

  const FastaSeqPtr& getFastaSeqPtr(int prot_idx) const {
    return seqs_[prot_idx];
  }

  // the number of residues indexed
  size_t getResidueNum() const;

  int getProtLen(int prot_idx) const {
    return static_cast<int>(seq_start_[prot_idx + 1] - seq_start_[prot_idx]) -
           1;
  }

  // the mass of the residues before position pos of protein prot_idx (pos
  // from 0 to the protein's length, the latter giving the residue mass sum)
  double getPrefixMass(int prot_idx, int pos) const {
    return prefix_masses_[seq_start_[prot_idx] + pos];
  }

  // the letter a residue letter is indexed as (I as L)
  static char indexLetter(char c) { return c == 'I' ? 'L' : c; }

 private:
  // 'A'..'Z' -> 0..25 (I as L), any other character -> 26
  static constexpr uint32_t ALPHABET_SIZE = 27;

  static uint32_t charCode(char c);

  void buildIndex();

  FastaSeqPtrVec seqs_;
  // the letters of all proteins, a separator after each; seq_start_[i] is
  // where protein i begins
  std::string text_;
  std::vector<uint32_t> seq_start_;
  // prefix_masses_[seq_start_[i] + k]: the mass of the first k residues of
  // protein i (k = its length, the separator's slot, the residue mass sum)
  std::vector<double> prefix_masses_;
  // positions_[bucket_start_[c] .. bucket_start_[c + 1]) are the text
  // positions, increasing, where the K-mer with code c starts
  std::vector<uint32_t> bucket_start_;
  std::vector<uint32_t> positions_;
};

using TagIndexPtr = std::shared_ptr<TagIndex>;

}  // namespace toppic

#endif
