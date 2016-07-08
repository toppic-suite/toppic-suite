#ifndef PROT_DECONV_ARGUMENT_HPP_
#define PROT_DECONV_ARGUMENT_HPP_

#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>

#include "boost/program_options.hpp" 

#include "base/logger.hpp"
#include "base/xml_dom_document.hpp"
#include "base/string_util.hpp"

namespace prot {

class DeconvArgument {
 public:
  DeconvArgument();
  static void outputArguments(std::ostream &output, 
                              std::map<std::string, std::string> arguments);
  bool parse(int argc, char* argv[]);
  std::map<std::string,std::string> getArguments(){return arguments_;}

 private:
  void initArguments();
  void setArgumentsByConfigFile(const std::string &file_name);
  bool validateArguments();
  void showUsage(boost::program_options::options_description &desc);
  std::map<std::string,std::string> arguments_;
};

}  // namespace prot

#endif /* ARGUMENT_HPP_ */
