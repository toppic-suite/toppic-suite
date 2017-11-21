//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include "protein_db.hpp"

namespace prot {
namespace suffix {

void ProteinDatabase::addProteinID(const std::string & proteinName) {
  std::string res = proteinName;
  int index = proteinName.find_first_of(" ", 0);
  if (index != -1) res = proteinName.substr(0, index);
  proteinID.push_back(res.substr(1));
}

}  // namespace suffix
}  // namespace prot
