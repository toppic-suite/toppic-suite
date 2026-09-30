// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#ifndef TOPPIC_CONSOLE_TOPLIB_BUILD_ARGUMENT_HPP_
#define TOPPIC_CONSOLE_TOPLIB_BUILD_ARGUMENT_HPP_

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "boost/program_options.hpp"

namespace toppic {

class ToplibBuildArgument {
 public:
  ToplibBuildArgument();

  static std::map<std::string, std::string> initArguments();

  static void outputArguments(std::ostream& output,
                              std::map<std::string, std::string> arguments);

  // false when the program should stop: an error, or the help was printed
  // (helpShown())
  bool parse(int argc, char* argv[]);

  bool helpShown() const { return help_shown_; }

  std::map<std::string, std::string> getArguments() { return arguments_; }

 private:
  bool validateArguments();

  void showUsage(boost::program_options::options_description& desc);

  std::map<std::string, std::string> arguments_;

  bool help_shown_ = false;
};

}  // namespace toppic

#endif
