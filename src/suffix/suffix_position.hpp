//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_SUFFIX_POSITION_HPP
#define PROT_SUFFIX_POSITION_HPP

#include <vector>
#include <string>

#include "protein_db.hpp"

namespace prot {
namespace suffix {

class SuffixPosition {
 public:
  SuffixPosition(int seq, int pos):
      seqNum(seq),
      posInSeq(pos),
      peptideStartPos(0),
      peptideEndPos(0) {}

  SuffixPosition(int seq, int pos, int pepStart):
      seqNum(seq),
      posInSeq(pos),
      peptideStartPos(pepStart),
      peptideEndPos(0) {}

  int getSeqNum() {return seqNum;}

  int getPosInSeq() {return posInSeq;}

  void setSuffixPos(int seq, int pos) {
    seqNum = seq;
    posInSeq = pos;
  }

  void setSeqNum(int seq) {seqNum = seq;}

  void setPosInSeq(int pos) {posInSeq = pos;}

  void setPeptideStartPos(int pepStart) {peptideStartPos = pepStart;}

  int getPeptideStartPos() {return peptideStartPos;}

  int getPeptideEndPos() {return peptideEndPos;}

  std::string toString();

 private:
  int seqNum;  // seq index starts from 0
  int posInSeq;
  int peptideStartPos;  // index starts from 0
  int peptideEndPos;  // index starts from 0, peptide ranges from [peptideStartPos, peptideEndPos)
};
}  // namespace suffix
}  // namespace prot
#endif
