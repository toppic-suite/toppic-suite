#include <iostream>
#include <sstream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "base/logger.hpp"
#include "base/file_util.hpp"

#ifdef PROT_LINUX

#include <unistd.h>

#endif

namespace fs = boost::filesystem;

namespace prot {

#ifdef PROT_LINUX
//linux
  std::string getExecutiveDir(const std::string &argv_0) {
    int buffer_size = 1024;
    char* buffer = new char[buffer_size];
    size_t len = readlink ("/proc/self/exe", buffer, buffer_size);
    std::string file_name;
    if (len >= 0) {
      buffer[len] = '\0';
      file_name = std::string (buffer);
    } 
    else {
      LOG_ERROR("Can not find executive directory!");
      return file_name = "";
    }
    delete buffer;
    fs::path full_path(file_name);
    std::string exe_dir = full_path.remove_filename().string();
    return exe_dir;
  }

#else 

//windows
  std::string getExecutiveDir(const std::string &argv_0) {
    fs::path full_path;
    full_path = fs::system_complete( fs::path(argv_0));

    std::string exe_dir = full_path.remove_filename().string();
    return exe_dir;
  }

#endif

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
  //std::cout << "parent dir " << parent_dir << std::endl;
  return parent_dir;
}

void createFolder(const std::string &folder_name){
  //std::cout<<"create folder "<<folder_name<<std::endl;
  boost::filesystem::path path(folder_name);
  boost::filesystem::create_directories(path);
}

void copyFile(const std::string &file_name, 
              const std::string &to_file, bool over_write){
  //std::cout<<"copy file: "<<file_name<<" to "<<to_file<<std::endl;
  boost::filesystem::path from_path(file_name);
  boost::filesystem::path to_path(to_file);
  if(!boost::filesystem::exists(from_path)){
    LOG_ERROR( "source file "<<file_name<<" does not exist!");
    return;
  }
  if(boost::filesystem::exists(to_path)){
    if(over_write){
      boost::filesystem::remove(to_path);
    }
    else{
      return;
    }
  }

  boost::filesystem::copy(from_path,to_path);
}

}
