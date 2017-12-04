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


#ifndef PROT_BASE_DB_BLOCK_HPP_
#define PROT_BASE_DB_BLOCK_HPP_

#include <memory>
#include <vector>
#include <string>

namespace prot {

class DbBlock;
typedef std::shared_ptr<DbBlock> DbBlockPtr;
typedef std::vector<DbBlockPtr> DbBlockPtrVec;

class DbBlock {
 public:
  DbBlock(int block_index, int seq_index):
      block_index_(block_index),
      seq_index_(seq_index) {}

  int getBlockIdx() {return block_index_;}

  int getSeqIdx() {return seq_index_;}

  static DbBlockPtrVec readDbBlockIndex(const std::string &index_file_name);

 private:
  int block_index_;
  int seq_index_;
};

}  // namespace prot
#endif
