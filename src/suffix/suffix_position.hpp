//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef PROT_SUFFIX_SUFFIX_POSITION_HPP
#define PROT_SUFFIX_SUFFIX_POSITION_HPP

#include <vector>
#include <string>
#include <memory>

#include "protein_db.hpp"

namespace prot {

namespace suffix {

class SuffixPosition {
 public:
  SuffixPosition(int seq, int pos):
      seq_num_(seq),
      pos_in_seq_(pos),
      peptide_start_pos_(0),
      peptide_end_pos_(0) {}

  SuffixPosition(int seq, int pos, int pep_start):
      seq_num_(seq),
      pos_in_seq_(pos),
      peptide_start_pos_(pep_start),
      peptide_end_pos_(0) {}

  int getSeqNum() {return seq_num_;}

  int getPosInSeq() {return pos_in_seq_;}

  void setSuffixPos(int seq, int pos) {
    seq_num_ = seq;
    pos_in_seq_ = pos;
  }

  void setSeqNum(int seq) {seq_num_ = seq;}

  void setPosInSeq(int pos) {pos_in_seq_ = pos;}

  void setPeptideStartPos(int pep_start) {peptide_start_pos_ = pep_start;}

  int getPeptideStartPos() {return peptide_start_pos_;}

  int getPeptideEndPos() {return peptide_end_pos_;}

 private:
  int seq_num_;  // seq index starts from 0

  int pos_in_seq_;

  int peptide_start_pos_;  // index starts from 0

  int peptide_end_pos_;  // index starts from 0, peptide ranges from [peptideStartPos, peptideEndPos)
};

typedef std::shared_ptr<SuffixPosition> SuffixPosPtr;

}  // namespace suffix

}  // namespace prot

#endif
