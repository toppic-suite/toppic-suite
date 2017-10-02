#define BOOST_TEST_DYN_LINK

#include <iostream>
#include <fstream>

#include <catch.hpp>

#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"

using namespace prot;

TEST_CASE("fasta reader") {
  std::ofstream test_file;
  test_file.open("test.fa");
  test_file << ">sp|test|test test_desc" << std::endl;
  test_file << "MSGRGKBGGXKGJLGAKG" << std::endl;
  test_file.close();

  BaseData::init("./");
  FastaReader f_reader("test.fa");
  FastaSeqPtr seq = f_reader.getNextSeq();
  REQUIRE(seq->getName() == "sp|test|test");
  REQUIRE(seq->getString(seq->getAcidPtmPairVec()) == "MSGRGKDGGAKGILGAKG");
  f_reader.close();
}

