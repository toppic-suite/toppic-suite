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


#ifndef PROT_BASE_FILE_UTIL_HPP_
#define PROT_BASE_FILE_UTIL_HPP_

#include <string>

#ifndef BOOST_SYSTEM_NO_DEPRECATED
#define BOOST_SYSTEM_NO_DEPRECATED 1
#endif

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

namespace prot {

namespace file_util {

std::string getFileSeparator();

std::string getExecutiveDir(const std::string &argv_0);

inline std::string getResourceDirName() {return "toppic_resources";}

std::string basename(const std::string &s);

std::string directory(const std::string &s);

void createFolder(const std::string &folder_name);

void copyFile(const std::string &file_name, const std::string &path,
              bool over_write);

bool copyDir(boost::filesystem::path const & source,
             boost::filesystem::path const & destination);

void delDir(const std::string &path);

void delFile(const std::string &path);

void cleanTopmgDir(const std::string &fa_path, const std::string & sp_path);

void cleanToppicDir(const std::string &fa_path, const std::string & sp_path);

void cleanTempFiles(const std::string & sp_path, 
                    const std::string & ext_prefix);

}  // namespace file_util

}  // namespace prot
#endif
