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


#ifndef PROT_PRSM_SIMPLE_PRSM_STR_HPP_
#define PROT_PRSM_SIMPLE_PRSM_STR_HPP_

#include <memory>
#include <vector>
#include <string>

namespace prot {

class SimplePrsmStr;

typedef std::shared_ptr<SimplePrsmStr> SimplePrsmStrPtr;

class SimplePrsmStr {
 public:
  explicit SimplePrsmStr(const std::vector<std::string> &str_vec);

  std::vector<std::string> getStrVec() {return str_vec_;}

  std::string getFileName() {return file_name_;}

  int getSpectrumId() {return spectrum_id_;}

  std::string getSeqName() {return seq_name_;}

  std::string getSeqDesc() {return seq_desc_;}

  double getScore() {return score_;}

  static bool cmpScoreDec(const SimplePrsmStrPtr &a, const SimplePrsmStrPtr &b) {
    if (a->getScore() == b->getScore()) {
      return a->getSeqName() < b->getSeqName();
    } else {
      return a->getScore() > b->getScore();
    }
  }

 private:
  std::vector<std::string> str_vec_;

  std::string file_name_;

  int spectrum_id_;

  double score_;

  std::string seq_name_;

  std::string seq_desc_;
};

typedef std::vector<SimplePrsmStrPtr> SimplePrsmStrPtrVec;

}  // namespace prot

#endif

