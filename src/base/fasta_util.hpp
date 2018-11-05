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


#ifndef PROT_BASE_FASTA_UTIL_HPP_
#define PROT_BASE_FASTA_UTIL_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#ifndef BOOST_SYSTEM_NO_DEPRECATED
#define BOOST_SYSTEM_NO_DEPRECATED 1
#endif

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

#include "htslib/faidx.h"
#include "base/fasta_reader.hpp"

namespace prot {

namespace fasta_util {

void generateShuffleDb(const std::string &file_name,
                       const std::string &target_decoy_file_name);

void dbSimplePreprocess(const std::string &ori_db_file_name,
                        const std::string &db_file_name);

void dbPreprocess(const std::string &ori_db_file_name,
                  const std::string &db_file_name,
                  bool decoy, int block_size);

int countProteinNum(const std::string &fasta_file);

}  // namespace fasta_util

}  // namespace prot

#endif
