//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#include <string>

#include "common/util/logger.hpp"
#include "seq/fasta_reader.hpp"

namespace toppic {

FastaReader::FastaReader(const std::string &file_name) {
  file_name_ = file_name;
  input_.open(file_name.c_str(), std::ios::in);
  if (!input_.is_open()) {
    LOG_ERROR("Fasta file  " << file_name << " does not exist.");
    exit(EXIT_FAILURE);
  }
  // if file is empty
  if (input_.peek() == std::ifstream::traits_type::eof()) {
    input_.close();
    return;
  }
  std::getline(input_, ori_name_);

  while (ori_name_.length() <= 1) {
    std::getline(input_, ori_name_);
  }
}

FastaReader::~FastaReader() {
  if (input_.is_open()) {
    input_.close();
  }
}

FastaSeqPtr FastaReader::getNextSeq() {
  if (!input_.is_open()) {
    return FastaSeqPtr(nullptr);
  }
  // get the letters of sequence
  std::string ori_seq;
  str_util::trim(ori_name_);
  std::string prot_name = ori_name_.substr(1, ori_name_.size() - 1);
  std::string line;
  while (std::getline(input_, line)) {
    if (line.length() >= 1 && line.substr(0, 1) == ">") {
      str_util::trim(line);
      ori_name_ = line;
      return std::make_shared<FastaSeq>(prot_name, ori_seq);
    }
    str_util::trim(line);
    ori_seq = ori_seq + line;
    if (ori_seq.size() >= 1000000) {
      LOG_ERROR("Protein sequences in " << file_name_ 
                << " are too long! Incorrect fasta file format!");
      exit(EXIT_FAILURE);
    }
  }
  input_.close();
  return std::make_shared<FastaSeq>(prot_name, ori_seq);
}

void FastaReader::close() {
  if (input_.is_open()) {
    input_.close();
  }
}

}  // namespace toppic
