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


#include "base/logger.hpp"
#include "base/fasta_reader.hpp"

namespace prot {

FastaReader::FastaReader(const std::string &file_name) {
  input_.open(file_name.c_str(), std::ios::in);
  if (!input_.is_open()) {
    LOG_ERROR( "fasta file  " << file_name << " does not exist.");
    throw "fasta file does not exist.";
  }
  std::getline(input_, ori_name_);

  while(ori_name_.length() <= 1){
    std::getline(input_, ori_name_);
  }  

}

FastaSeqPtr FastaReader::getNextSeq() {
  if (!input_.is_open()) {
    return FastaSeqPtr(nullptr);
  }

  // get the letters of sequence 
  std::string ori_seq;
  ori_name_ = StringUtil::trim(ori_name_);
  std::string prot_name = ori_name_.substr(1, ori_name_.size() - 1);
  std::string line;
  while (std::getline(input_, line)) {
    if (line.length() >= 1 && line.substr(0, 1) == ">") {
      ori_name_ = StringUtil::trim(line);
      return FastaSeqPtr(new FastaSeq(prot_name, ori_seq));
    }
    line = StringUtil::trim(line);
    ori_seq = ori_seq + line;
    if (ori_seq.size() >= 1000000) {
      LOG_ERROR("Protein sequences are too long! Incorrect fasta file!");
      throw("fasta file error");
    }
  }
  input_.close();
  return FastaSeqPtr(new FastaSeq(prot_name, ori_seq));
}

void FastaReader::close() {
  input_.close();
}

}
