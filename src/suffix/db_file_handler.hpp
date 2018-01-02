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


#ifndef PROT_SUFFIX_DB_FILEHANDLER_HPP
#define PROT_SUFFIX_DB_FILEHANDLER_HPP

#include <fstream>
#include <iostream>
#include <string>

#include "protein_db.hpp"

namespace prot {

namespace suffix {

class DatabaseFileHandler {
 public:
  ProteinDBPtr loadDatabase(const std::string & db_file);

 private:
  std::string handleUndefinedCharacter(std::string text);
};

}  // namespace suffix

}  // namespace prot
#endif
