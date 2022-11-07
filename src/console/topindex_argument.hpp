//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_TOPINDEX_ARGUMENT_HPP_
#define TOPPIC_TOPINDEX_ARGUMENT_HPP_

#include <map>
#include <iostream>
#include <memory>
#include <vector>
#include <string>

#include "boost/program_options.hpp" 

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"

namespace toppic {

class TopIndexArgument {
 public:
  TopIndexArgument();

  static std::map<std::string, std::string> initArguments();

  static void outputArguments(std::ostream &output, 
		  	      const std::string &sep,
                              std::map<std::string, std::string> arguments);

  bool parse(int argc, char* argv[]);

  std::map<std::string, std::string> getArguments() {return arguments_;}

 private:

  bool validateArguments();

  void showUsage(boost::program_options::options_description &desc);

  std::map<std::string, std::string> arguments_;
};

}  // namespace toppic

#endif /* ARGUMENT_HPP_ */
