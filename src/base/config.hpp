/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_CONF_HPP_
#define PROT_BASE_CONF_HPP_

#include <string>

namespace prot {

class Config {
  static std::string getConfigDir() {return "toppic_resources/conf";}

  static std::string getConfigFileName() {return "configuration.xml";}

  static void initConf(const std::string &exe_dir);
};

}

#endif

