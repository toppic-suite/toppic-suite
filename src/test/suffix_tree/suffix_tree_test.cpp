// Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <iostream>
#include <fstream>

#include <catch.hpp>

#include "suffix/db_file_handler.hpp"
#include "suffix/suffix_tree.hpp"


using namespace toppic::suffix;

TEST_CASE("suffix tree") {
  std::ofstream test_file;
  test_file.open("test.fa");
  test_file << ">sp|test1|test1 test1_desc" << std::endl;
  test_file << "MSGRGKBGGXKGJLGAKG" << std::endl;
  test_file << ">sp|test2|test2 test2_desc" << std::endl;
  test_file << "MSKGRGKGKLRXKGJLGAKG" << std::endl;

  test_file.close();

  std::shared_ptr<DatabaseFileHandler> df = std::make_shared<DatabaseFileHandler>();
  ProteinDBPtr pd = df->loadDatabase("test.fa");
  REQUIRE(pd->getsize() == 2);
  REQUIRE(pd->getProteinID(0) == "sp|test1|test1");
  REQUIRE(pd->getProteinDesc(1) == "test2_desc");
  REQUIRE(pd->getSequence() == "MSGRGKAGGAKGAIGAKG#MSKGRGKGKIRAKGAIGAKG$");

  std::shared_ptr<SuffixTree> st = std::make_shared<SuffixTree>(pd->getSequence(), pd);
  st->init();

  std::vector<SuffixPosPtr> startPosList = st->search("GRG");

  REQUIRE(startPosList.size() == 2);
  REQUIRE(startPosList[0]->getSeqNum() == 0);
  REQUIRE(startPosList[0]->getPosInSeq() == 2);
  REQUIRE(startPosList[1]->getSeqNum() == 1);
  REQUIRE(startPosList[1]->getPosInSeq() == 3);
}

