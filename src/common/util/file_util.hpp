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


#ifndef TOPPIC_COMMON_UTIL_FILE_UTIL_HPP_
#define TOPPIC_COMMON_UTIL_FILE_UTIL_HPP_

#include <string>

namespace toppic {

namespace file_util {

std::string getExecutiveDir(const std::string &argv_0);

std::string getFileSeparator();

std::string basename(const std::string &s);

std::string directory(const std::string &s);

std::string absoluteDir(const std::string &s);

std::string absoluteName(const std::string &s);

void createFolder(const std::string &folder_name);

void copyFile(const std::string &file_name, const std::string &path,
              bool over_write);

bool copyDir(const std::string &source, 
             const std::string &destination);

bool exists(const std::string &path);

void delDir(const std::string &path);

void delFile(const std::string &path);

void rename(const std::string &ori_name, 
            const std::string &new_name);

void cleanPrefix(const std::string & ref_name, 
                 const std::string & prefix);

void cleanTempFiles(const std::string & ref_name, 
                    const std::string & ext_prefix);

inline std::string getToppicResourceDirName() {return "toppic_resources";}

inline std::string getEtcDirName() {return "../../etc/toppic";}

std::string getResourceDir(const std::string &exec_dir);

}  // namespace file_util

}  // namespace toppic
#endif
