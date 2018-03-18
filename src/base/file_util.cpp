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


#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "base/logger.hpp"
#include "base/file_util.hpp"

namespace fs = boost::filesystem;

namespace prot {

namespace file_util {

std::string getFileSeparator() {
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  return "\\";
#else
  return "/";
#endif
}

std::string getExecutiveDir(const std::string &argv_0) {
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  LPSTR lpFilePart;
  char file_name[MAX_PATH];
  SearchPath(NULL, argv_0.c_str(), ".exe", MAX_PATH, file_name, &lpFilePart);
#else
  int buffer_size = 1024;
  char* buffer = new char[buffer_size];
  size_t len = readlink("/proc/self/exe", buffer, buffer_size);
  std::string file_name;
  buffer[len] = '\0';
  file_name = std::string(buffer);
  delete[] buffer;
#endif
  fs::path full_path(file_name);
  std::string exe_dir = full_path.remove_filename().string();
  return exe_dir;
}

std::string basename(const std::string &s) {
  size_t dot_pos = s.find_last_of(".");
  if (dot_pos < s.length()) {
    return s.substr(0, dot_pos);
  }
  return s;
}

std::string directory(const std::string &s) {
  fs::path path(s);
  std::string parent_dir = path.parent_path().string();
  return parent_dir;
}

void createFolder(const std::string &folder_name) {
  boost::filesystem::path path(folder_name);
  boost::filesystem::create_directories(path);
}

void copyFile(const std::string &file_name,
                        const std::string &to_file, bool over_write) {
  boost::filesystem::path from_path(file_name);
  boost::filesystem::path to_path(to_file);
  if (!boost::filesystem::exists(from_path)) {
    LOG_ERROR("source file " << file_name << " does not exist!");
    return;
  }

  if (boost::filesystem::exists(to_path)) {
    if (over_write) {
      boost::filesystem::remove(to_path);
    } else {
      return;
    }
  }

  boost::filesystem::copy(from_path, to_path);
}

bool copyDir(boost::filesystem::path const & source,
                       boost::filesystem::path const & destination) {
  try {
    if (!fs::exists(source) || !fs::is_directory(source)) {
      return false;
    }
    if (fs::exists(destination)) {
      return false;
    }

    if (!fs::create_directory(destination)) {
      std::cerr << "Unable to create destination directory"
          << destination.string() << std::endl;
      return false;
    }
  }
  catch(fs::filesystem_error const & e) {
    std::cerr << e.what() << std::endl;
    return false;
  }

  for (fs::directory_iterator file(source); file != fs::directory_iterator(); ++file) {
    try {
      fs::path current(file->path());
      if (fs::is_directory(current)) {
        if (!copyDir(current, destination / current.filename())) {
          return false;
        }
      } else {
        fs::copy_file(current, destination / current.filename());
      }
    }
    catch(fs::filesystem_error const & e) {
      std:: cerr << e.what() << std::endl;
    }
  }
  return true;
}

void delDir(const std::string &path) {
  fs::path dir(path);

  if (fs::exists(dir))
    fs::remove_all(dir);
}

void delFile(const std::string &path) {
  fs::path dir(path);

  if (fs::exists(dir))
    fs::remove(dir);
}

void clean_prefix(const fs::path & sp, const std::string & prefix) {
  fs::directory_iterator end_iter;
  for (fs::directory_iterator dir_iter(absolute(sp).parent_path());
       dir_iter != end_iter ; ++dir_iter) {
    std::string filename = dir_iter->path().string();
    std::replace(filename.begin(), filename.end(), '\\', '/');
    if (filename.compare(0, prefix.length(), prefix) == 0) {
      if (!fs::is_directory(fs::status(dir_iter->path())))
        fs::remove(dir_iter->path());
    }
  }
}

void cleanDir(const std::string &fa_path, const std::string & sp_path) {
  fs::path fa(fa_path);
  fs::path sp(sp_path);
  std::string fa_base = absolute(fa).string();
  std::replace(fa_base.begin(), fa_base.end(), '\\', '/');
  std::string sp_base = basename(absolute(sp).string());
  std::replace(sp_base.begin(), sp_base.end(), '\\', '/');

  clean_prefix(fa, fa_base + "_");
  clean_prefix(sp, sp_base + ".msalign_");
  delFile(absolute(sp).string() + "_index");
  delFile(sp_base + ".ZERO_PTM");
  clean_prefix(sp, sp_base + ".ZERO_PTM_");
  clean_prefix(sp, sp_base + ".ZERO_FILTER_");
  delFile(sp_base + ".ONE_PTM");
  clean_prefix(sp, sp_base + ".ONE_PTM_");
  delFile(sp_base + ".DIAG_FILTER");
  clean_prefix(sp, sp_base + ".DIAG_FILTER_");
  delFile(sp_base + ".PTM");
  clean_prefix(sp, sp_base + ".PTM_");
  delFile(sp_base + ".TOP");
  delFile(sp_base + ".TOP_PRE");
  delFile(sp_base + ".GRAPH_FILTER");
  clean_prefix(sp, sp_base + ".GRAPH_ALIGN_");
  clean_prefix(sp, sp_base + ".GRAPH_FILTER_");
  clean_prefix(sp, sp_base + ".GRAPH");
  clean_prefix(sp, sp_base + ".VAR1_");
  clean_prefix(sp, sp_base + ".VAR2_");
  delFile(sp_base + ".GRAPH_ALIGN");
  delFile(sp_base + ".CUTOFF_RESULT_SPEC");
  delFile(sp_base + ".CUTOFF_RESULT_FORM");
  delFile(sp_base + ".LOCAL_RESULT_SPEC");
  delFile(sp_base + ".LOCAL_RESULT_FORM");
  clean_prefix(sp, sp_base + ".EVALUE_");
  clean_prefix(sp, sp_base + ".PROTEOFORM_");
  clean_prefix(sp, sp_base + ".PROT_");
  delFile(sp_base + ".EVALUE");
  delFile(sp_base + ".RAW_RESULT");
  delFile(sp_base + ".CLUSTERS");
  delFile(sp_base + ".FORM_RESULT");
  delFile(sp_base + ".FORM_FILTER_RESULT");
}

}  // namespace file_util

}  // namespace prot
