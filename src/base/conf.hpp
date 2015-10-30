/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_CONF_HPP_
#define PROT_BASE_CONF_HPP_

#include <string>

#define CONFIG_DIR  "toppic_resources/conf"
#define CONFIG_FILE_NAME "configuration.xml"

namespace prot {


void initConf(const std::string &exe_dir);

}
#endif

