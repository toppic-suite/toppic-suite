/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_BASE_DATA_HPP_
#define PROT_BASE_BASE_DATA_HPP_

#include <string>

namespace prot {

class BaseData {
  static std::string getBaseDataDir() {return "toppic_resources/base_data";}

  static std::string getBaseDataConfigFileName() {return "base_data_config.xml";}

  static void init(const std::string &exe_dir);
};

}

#endif

