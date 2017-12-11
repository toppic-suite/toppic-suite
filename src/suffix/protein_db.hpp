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

#ifndef PROT_SUFFIX_DB_HPP
#define PROT_SUFFIX_DB_HPP

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <vector>

namespace prot {

namespace suffix {

class ProteinDatabase {
 public:
  ProteinDatabase(): sequence("") {}

  size_t getsize() {return seqSet.size();}

  size_t getSeqLength(int seqIndex) {return seqLen[seqIndex];}

  std::string getIndividualSeq(int index) {return seqSet[index];}

  void addIndividualSeq(const std::string & seq) {
    seqLen.push_back(seq.length());
    seqSet.push_back(seq);
  }

  std::string getSequence() {return sequence;}

  void setSequence(const std::string & seq) {sequence = seq;}

  std::string getProteinID(int index) {return proteinID[index];}

  std::string getProteinDesc(int index) {return proteinDesc[index];}

  void addProteinID(const std::string & proteinName);

 private:
  std::string sequence;  // combined sequence

  std::vector<std::string> proteinID;

  std::vector<std::string> proteinDesc;

  std::vector<std::string> seqSet;  // store each individual sequence

  std::vector<size_t> seqLen;
};

}  // namespace suffix

}  // namespace prot
#endif
