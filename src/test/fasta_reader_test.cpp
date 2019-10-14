//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include <iostream>
#include <fstream>

#include <catch.hpp>

#include "common/base/base_data.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/fasta_sub_util.hpp"

using namespace toppic;

TEST_CASE("fasta reader") {
  std::ofstream test_file;
  test_file.open("test.fa");
  test_file << ">sp|test|test test_desc" << std::endl;
  test_file << "MSGRGKBGGXKGJLGAKG" << std::endl;
  test_file.close();

  base_data::init();
  FastaReader f_reader("test.fa");
  FastaSeqPtr seq = f_reader.getNextSeq();
  f_reader.close();
  REQUIRE(seq->getName() == "sp|test|test");
  REQUIRE(seq->getRawSeq() == "MSGRGKBGGXKGJLGAKG");

  // test break_seq
  FastaSubSeqPtrVec sub_seqs = fasta_sub_util::breakSeq(seq, 10);
  REQUIRE(sub_seqs[0]->getName() == "sp|test|test");
  REQUIRE(sub_seqs[0]->getRawSeq() == "MSGRGKDGGA");
  REQUIRE(sub_seqs[1]->getRawSeq() == "KDGGAKGILG");
  REQUIRE(sub_seqs[2]->getRawSeq() == "KGILGAKG");
}

