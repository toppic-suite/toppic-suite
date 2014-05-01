#ifndef PROT_FILE_UTIL_HPP_
#define PROT_FILE_UTIL_HPP_

#include <string>

//#define LINUX

namespace prot {

#define FILE_SEPARATOR "/"

std::string getExecutiveDir(const std::string &argv_0);

std::string basename(const std::string &s);

std::string directory(const std::string &s);

void createFolder(const std::string &folder_name);

void copyFile(const std::string &file_name, const std::string &path, 
              bool over_write);

}

#endif
