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

#include <algorithm>
#include <string>

#include "db_file_handler.hpp"

namespace toppic {

namespace suffix {

ProteinDBPtr DatabaseFileHandler::loadDatabase(const std::string & db_file) {
  ProteinDBPtr database = std::make_shared<ProteinDatabase>();
  std::ifstream fis;
  fis.open(db_file.c_str());
  if (!fis) {
    std::cerr << "Unable to open Database " << db_file << std::endl;
    exit(EXIT_FAILURE);
  }

  bool isFirstSeq = true;
  std::string sequence;
  std::string individualSeq;
  std::string line;

  while (std::getline(fis, line)) {
    if (line[0] == '>') {
      database->addProteinID(line);
      if (!isFirstSeq) {
        sequence.append("#");
        database->addIndividualSeq(individualSeq);
        individualSeq = "";
      }
      isFirstSeq = false;
    } else {
      line = handleUndefinedCharacter(line);
      individualSeq.append(line);
      sequence.append(line);
    }
  }
  sequence.append("$");
  sequence.erase(std::remove(sequence.begin(), sequence.end(), '\r'), sequence.end());
  std::replace(sequence.begin(), sequence.end(), 'L', 'I');
  database->setSequence(sequence);
  database->addIndividualSeq(individualSeq);

  sequence = "";
  individualSeq = "";
  fis.close();
  return database;
}

std::string DatabaseFileHandler::handleUndefinedCharacter(std::string text) {
  std::replace(text.begin(), text.end(), 'B', 'A');
  std::replace(text.begin(), text.end(), 'J', 'A');
  std::replace(text.begin(), text.end(), 'O', 'A');
  std::replace(text.begin(), text.end(), 'U', 'A');
  std::replace(text.begin(), text.end(), 'X', 'A');
  std::replace(text.begin(), text.end(), 'Z', 'A');
  return text;
}

}  // namespace suffix

}  // namespace toppic
