// Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_CONSOLE_TOPDIA_ARGUMENT_HPP_
#define TOPPIC_CONSOLE_TOPDIA_ARGUMENT_HPP_

#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "topfd/common/topfd_para.hpp"
#include "topdia/common/topdia_para.hpp"

namespace toppic {

class Argument {
 public:
  Argument();

  bool parse(int argc, char *argv[]);

  static TopfdParaPtr getTopfdParaPtrForTopdia();

  TopfdParaPtr getTopfdParaPtr() { return topfd_para_ptr_; }

  TopdiaParaPtr getTopdiaParaPtr() { return topdia_para_ptr_; }

  std::vector<std::string> getSpecFileList() { return spec_file_list_; };

 private:
  void initArguments();

  void setArgumentsByConfigFile(const std::string &file_name);

  bool validateArguments();

  void showUsage(boost::program_options::options_description &desc);

  TopfdParaPtr topfd_para_ptr_;

  TopdiaParaPtr topdia_para_ptr_;

  std::vector<std::string> spec_file_list_;
};

}  // namespace toppic

#endif