#ifndef PROT_BASE_FASTA_UTIL_HPP_
#define PROT_BASE_FASTA_UTIL_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

#include "htslib/faidx.h"
#include "base/fasta_reader.hpp"

namespace prot {

class FastaUtil {
 public:
  static void generateShuffleDb(const std::string &file_name, 
                                const std::string &target_decoy_file_name);

  static void dbPreprocess(const std::string &ori_db_file_name, 
                           const std::string &db_file_name, 
                           bool decoy, int block_size);

  static int countProteinNum(const std::string &fasta_file);
};

}  //namepace prot

#endif
