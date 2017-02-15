// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "base/db_block.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "base/fasta_reader.hpp"

namespace prot {

DbBlock::DbBlock(int block_index, int seq_index) {
  block_index_ = block_index;
  seq_index_ = seq_index;
}

DbBlockPtrVec DbBlock::readDbBlockIndex(const std::string &db_file_name) {
  DbBlockPtrVec db_block_ptr_vec;
  std::ifstream input;
  std::string index_file_name = db_file_name + "_block_index";
  input.open(index_file_name.c_str(), std::ios::in);
  if (!input.is_open()) {
    LOG_ERROR( "index file  " << index_file_name << " does not exist.");
    throw "index file does not exist.";
  }
  std::string line;
  std::vector<std::string> strs;
  while (std::getline(input, line)) {
    strs = StringUtil::split(line, '\t');
    int block_index = std::stoi(strs[0]);
    int seq_index = std::stoi(strs[1]);
    LOG_DEBUG("block " << block_index << " seq " << seq_index);
    DbBlockPtr ptr(new DbBlock(block_index, seq_index));
    db_block_ptr_vec.push_back(ptr);
  }
  input.close();
  return db_block_ptr_vec;
}

}
