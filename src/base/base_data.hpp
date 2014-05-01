/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_DATA_HPP_
#define PROT_BASE_DATA_HPP_

#include <string>

#define CONFIG_DIR  "conf"
#define CONFIG_FILE_NAME "configuration.xml"

namespace prot {


void initBaseData(const std::string &exe_dir);

}
#endif

