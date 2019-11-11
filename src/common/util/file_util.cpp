//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include <string>

#ifndef BOOST_SYSTEM_NO_DEPRECATED
#define BOOST_SYSTEM_NO_DEPRECATED 1
#endif

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"

namespace fs = boost::filesystem;

namespace toppic {

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

std::string getResourceDir(const std::string &exec_dir) {
  std::string resource_dir = exec_dir + getFileSeparator() + getToppicResourceDirName();
  if (exists(resource_dir)) {
    return resource_dir;
  }
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
#else
  std::string etc_dir = exec_dir + getFileSeparator() + getEtcDirName();
  if (exists(etc_dir)) {
    return etc_dir;
  }
  LOG_ERROR("The resource directory " << etc_dir << " does not exist!"); 
#endif
  LOG_ERROR("The resource directory " << resource_dir << " does not exist!"); 
  exit(EXIT_FAILURE);
}

std::string basename(const std::string &s) {
  size_t dot_pos = s.find_last_of(".");
  if (dot_pos < s.length()) {
    return s.substr(0, dot_pos);
  }
  return s;
}

std::string absoluteDir(const std::string &s) {
  fs::path path(s);
  std::string parent_dir = absolute(path).parent_path().string();
  return parent_dir;
}

std::string absoluteName(const std::string &s) {
  fs::path path(s);
  std::string absolute_name = absolute(path).string();
  return absolute_name;
}

std::string directory(const std::string &s) {
  fs::path path(s);
  std::string parent_dir = path.parent_path().string();
  return parent_dir;
}

void createFolder(const std::string &folder_name) {
  fs::path path(folder_name);
  fs::create_directories(path);
}

void createLink(const std::string &a_link, const std::string &a_dir, const std::string &b) {
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  copyDir(a_dir, b);
#else
  fs::path path_a(a_link);
  fs::path path_b(b);
  fs::create_symlink(path_a, path_b);
#endif
}

void copyFile(const std::string &file_name,
              const std::string &to_file, bool over_write) {
  fs::path from_path(file_name);
  fs::path to_path(to_file);
  if (!fs::exists(from_path)) {
    LOG_ERROR("The source file " << file_name << " does not exist!");
    return;
  }

  if (over_write) {
    fs::copy_file(from_path, to_path, fs::copy_option::overwrite_if_exists);
  }
  else {
    fs::copy_file(from_path, to_path, fs::copy_option::fail_if_exists);
  }
}

bool copyDir(const std::string &src_name,
             const std::string &des_name) {
  fs::path source(src_name);
  fs::path destination(des_name);
  try {
    if (!fs::exists(source) || !fs::is_directory(source)) {
      LOG_ERROR("The source folder " << source.string() << " does not exist!");
      return false;
    }
    if (fs::exists(destination)) {
      LOG_ERROR("The destination folder " << destination.string() 
          << " already exists. Fail to create the destination directory.");
      return false;
    }

    if (!fs::create_directory(destination)) {
      LOG_ERROR("Unable to create the destination folder "
                << destination.string());
      return false;
    }
  }
  catch(fs::filesystem_error const & e) {
    LOG_ERROR(e.what());
    return false;
  }

  for (fs::directory_iterator file(source); file != fs::directory_iterator(); ++file) {
    try {
      fs::path current(file->path());
      if (fs::is_directory(current)) {
        if (!copyDir(current.string(), (destination / current.filename()).string())) {
          return false;
        }
      } else {
        fs::copy_file(current, destination / current.filename());
      }
    }
    catch(fs::filesystem_error const & e) {
      LOG_ERROR(e.what());
    }
  }
  return true;
}

bool copyJsonDir(const std::string &src_name,
                 const std::string &des_name,
                 int id_base) {
  fs::path source(src_name);
  fs::path destination(des_name);
  try {
    if (!fs::exists(source) || !fs::is_directory(source)) {
      LOG_ERROR("The source folder " << source.string() << " does not exist!");
      return false;
    }
  }
  catch(fs::filesystem_error const & e) {
    LOG_ERROR(e.what());
    return false;
  }

  for (fs::directory_iterator file(source); file != fs::directory_iterator(); ++file) {
    try {
      fs::path current(file->path());
      std::string file_name = current.filename().string();
      std::string id_str = file_name.substr(8, file_name.length() - 3 - 8);
      //LOG_ERROR(file_name << " " << id_str);
      int new_id = std::stoi(id_str) + id_base;
      std::string new_name = "spectrum" + std::to_string(new_id) + ".js";
      fs::path des_file(des_name + getFileSeparator() + new_name);
      fs::copy_file(current, des_file);
    }
    catch(fs::filesystem_error const & e) {
      LOG_ERROR(e.what());
    }
  }
  return true;
}

bool exists(const std::string &path) {
  return fs::exists(path); 
}

void delDir(const std::string &path) {
  fs::path dir(path);
  if (fs::exists(dir))
    fs::remove_all(dir);
}

void delFile(const std::string &path) {
  fs::path file(path);
  if (fs::exists(file))
    fs::remove(file);
}

void rename(const std::string &ori_name, 
            const std::string &new_name) {
  fs::path ori_path(ori_name); 
  fs::path new_path(new_name);
  fs::rename(ori_path, new_path);

}

void moveFile(std::string &path_name, std::string &folder_name) {
  fs::path ori_path(path_name);
  std::string new_path_name = folder_name + getFileSeparator() + ori_path.filename().string();
  bool over_write = true;
  copyFile(path_name, new_path_name, over_write); 
  delFile(path_name);
}

void cleanPrefix(const std::string & ref_name, 
                 const std::string & prefix) {
  fs::path ref_path(ref_name);
  fs::path ref_dir = absolute(ref_path).parent_path(); 
  fs::directory_iterator end_iter;
  for (fs::directory_iterator dir_iter(ref_dir); 
       dir_iter != end_iter ; ++dir_iter) {
    std::string file_name = dir_iter->path().string();
    std::replace(file_name.begin(), file_name.end(), '\\', '/');
    if (file_name.compare(0, prefix.length(), prefix) == 0) {
      if (!fs::is_directory(fs::status(dir_iter->path())))
        fs::remove(dir_iter->path());
    }
  }
}

void cleanTempFiles(const std::string & ref_name, 
                    const std::string & ext_prefix) {

  fs::path ref_path(ref_name);
  std::string ref_base = basename(absolute(ref_path).string());
  std::replace(ref_base.begin(), ref_base.end(), '\\', '/');
  cleanPrefix(ref_name, ref_base + "." + ext_prefix);
}

}  // namespace file_util

}  // namespace toppic
