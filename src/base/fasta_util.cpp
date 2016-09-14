// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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
//#include "base/db_block.hpp"
#include "base/fasta_util.hpp"

namespace prot {

void FastaUtil::generateShuffleDb(const std::string &file_name, 
                                  const std::string &target_decoy_file_name) {
  std::ofstream output;
  output.open(target_decoy_file_name.c_str(), std::ios::out);
  FastaReader reader(file_name);

  FastaSeqPtr seq_info = reader.getNextSeq();
  std::mt19937 r{std::random_device{}()};
  r.seed(std::mt19937::default_seed);
  while (seq_info!=nullptr) {
    std::string name = seq_info->getName();
    std::string seq = seq_info->getRawSeq();
    StringPairVec str_pair_vec = seq_info->getAcidPtmPairVec();
    std::string desc = seq_info->getDesc();
    std::string decoy_name = "DECOY_" + name;
    std::string decoy_seq;
    /*
    if (seq.length() > 2) {
      std::string temp = seq.substr(2, seq.length() - 2);
      std::shuffle(temp.begin(), temp.end(), r);
      //std::random_shuffle(temp.begin(), temp.end());
      decoy_seq = seq.substr(0,2) + temp;
    }
    else {
      decoy_seq = seq;
    }
    */
    if (str_pair_vec.size() > 2) {
      std::shuffle(str_pair_vec.begin() + 2, str_pair_vec.end(), r);
      decoy_seq = FastaSeq::getString(str_pair_vec); 
    }
    else {
      decoy_seq = seq;
    }
    output << ">" << decoy_name << " " << desc <<  std::endl;
    output << decoy_seq << std::endl;
    output << ">" << name << " " << desc << std::endl;
    output << seq << std::endl;

    seq_info = reader.getNextSeq();
  }
  output.close();
}

void generateStandardDb(const std::string &ori_file_name, 
                        const std::string &st_file_name) {

  std::ifstream ori_db(ori_file_name);
  std::string line;
  std::ofstream standard_db;
  standard_db.open(st_file_name.c_str(), std::ios::out);

  while (std::getline(ori_db, line)) {
    if(line.length() > 0){
      standard_db << line << std::endl;
    }  
  }
  ori_db.close();
  standard_db.close();
}

void generateDbBlock(const std::string &db_file_name, 
                     int block_size) {
  int block_idx = 0;
  int seq_idx = 0;

  std::ofstream index_output;
  std::string index_file_name = db_file_name + "_block_index";
  index_output.open(index_file_name.c_str(), std::ios::out);
  std::ofstream block_output;
  std::string block_file_name = db_file_name + "_" + std::to_string(block_idx);
  block_output.open(block_file_name.c_str(), std::ios::out);

  FastaReader reader(db_file_name);
  FastaSeqPtr seq_info = reader.getNextSeq();
  index_output << block_idx << "\t" << seq_idx << std::endl;
  int seq_size = 0;
  while (seq_info!=nullptr) {
    std::string name = seq_info->getName();
    std::string desc = seq_info->getDesc();
    std::string seq = seq_info->getRawSeq();
    seq_size += seq.length();
    if (seq_size > block_size) {
      block_output.close();
      block_idx++;
      LOG_DEBUG("Database block " << block_idx << " size " << seq_size);
      index_output << block_idx << "\t" << seq_idx << std::endl;
      seq_size = 0;
      block_file_name = db_file_name + "_" + std::to_string(block_idx);
      block_output.open(block_file_name.c_str(), std::ios::out);
    }
    block_output << ">" << name << " " << desc << std::endl;
    block_output << seq << std::endl;
    seq_info = reader.getNextSeq();
    seq_idx++;
  }
  index_output.close();
  block_output.close();
}


void FastaUtil::dbPreprocess(const std::string &ori_db_file_name, 
                             const std::string &db_file_name, 
                             bool decoy, int block_size) {

  std::string standard_db_file_name = ori_db_file_name + "_standard";
  generateStandardDb(ori_db_file_name, standard_db_file_name);

  if (decoy) {
    generateShuffleDb(standard_db_file_name, db_file_name);
  }
  else {
    boost::filesystem::path ori_path(standard_db_file_name);
    boost::filesystem::path db_path(db_file_name);
    boost::filesystem::copy_file(ori_path, db_path, 
                                 boost::filesystem::copy_option::overwrite_if_exists);
  }
  generateDbBlock(db_file_name, block_size);
  fai_build(db_file_name.c_str());
}

int FastaUtil::countProteinNum(const std::string &fasta_file) {
  FastaReader reader(fasta_file);
  int cnt = 0;
  while (reader.getNextSeq() != nullptr) {
    cnt++;
  }
  reader.close();
  return cnt;
}

}

