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


#ifndef PROT_TOPMG_ARGUMENT_HPP_
#define PROT_TOPMG_ARGUMENT_HPP_

#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include <string>

#include <boost/program_options.hpp>

#include "base/logger.hpp"
#include "base/xml_dom_document.hpp"
#include "base/string_util.hpp"

namespace toppic {

class Argument {
 public:
  Argument();

  static void outputArguments(std::ostream &output,
                              std::map<std::string, std::string> arguments);

  static std::string outputCsvArguments(std::map<std::string, std::string> arguments);

  bool parse(int argc, char* argv[]);

  std::map<std::string, std::string> getArguments() {return arguments_;}

  std::vector<std::string> getSpecFileList() { return spec_file_list_;}
 private:
  void initArguments();

  bool validateArguments();

  void showUsage(boost::program_options::options_description &desc);

  std::map<std::string, std::string> arguments_;

  std::vector<std::string> spec_file_list_;
};

}  // namespace toppic

#endif /* ARGUMENT_HPP_ */
