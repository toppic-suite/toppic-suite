//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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
#include <memory>

namespace toppic {

namespace suffix {

class ProteinDatabase {
 public:
  ProteinDatabase(): seq_("") {}

  size_t getsize() {return seq_set_.size();}

  size_t getSeqLength(int seqIndex) {return seq_len_[seqIndex];}

  std::string getIndividualSeq(int index) {return seq_set_[index];}

  void addIndividualSeq(const std::string & seq) {
    seq_len_.push_back(seq.length());
    seq_set_.push_back(seq);
  }

  std::string getSequence() {return seq_;}

  void setSequence(const std::string & seq) {seq_ = seq;}

  std::string getProteinID(int index) {return protein_id_[index];}

  std::string getProteinDesc(int index) {return protein_desc_[index];}

  void addProteinID(const std::string & proteinName) {
    int index = proteinName.find_first_of(" ", 0);
    if (index != -1) {
      std::string res = proteinName.substr(0, index);
      protein_id_.push_back(res.substr(1));
      protein_desc_.push_back(proteinName.substr(index + 1));
    }
  }

 private:
  std::string seq_;  // combined sequence

  std::vector<std::string> protein_id_;

  std::vector<std::string> protein_desc_;

  std::vector<std::string> seq_set_;  // store each individual sequence

  std::vector<size_t> seq_len_;
};

typedef std::shared_ptr<ProteinDatabase> ProteinDBPtr;

}  // namespace suffix

}  // namespace toppic
#endif
