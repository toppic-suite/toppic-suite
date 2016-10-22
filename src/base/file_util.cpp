// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <iostream>
#include <sstream>

#include "base/logger.hpp"
#include "base/file_util.hpp"

#ifdef PROT_LINUX

#include <unistd.h>

#endif

namespace fs = boost::filesystem;

namespace prot {

#ifdef PROT_LINUX
  //linux
  std::string FileUtil::getFileSeparator() {
    return "/";
  }

  std::string FileUtil::getExecutiveDir(const std::string &argv_0){
    int buffer_size = 1024;
    char* buffer = new char[buffer_size];
    size_t len = readlink ("/proc/self/exe", buffer, buffer_size);
    std::string file_name;
    buffer[len] = '\0';
    file_name = std::string (buffer);
    delete[] buffer;
    fs::path full_path(file_name);
    std::string exe_dir = full_path.remove_filename().string();
    return exe_dir;
  }

#else
  //windows
  std::string FileUtil::getFileSeparator() {
    return "\\";
  }

  std::string FileUtil::getExecutiveDir(const std::string &argv_0){
    fs::path full_path;
    full_path = fs::system_complete( fs::path(argv_0));

    std::string exe_dir = full_path.remove_filename().string();
    return exe_dir;
  }
#endif

  std::string FileUtil::basename(const std::string &s){
    size_t dot_pos = s.find_last_of(".");
    if (dot_pos < s.length()){
      return s.substr(0, dot_pos);
    }
    return s;
  }

  std::string FileUtil::directory(const std::string &s){
    fs::path path(s);
    std::string parent_dir = path.parent_path().string();
    //std::cout << "parent dir " << parent_dir << std::endl;
    return parent_dir;
  }

  void FileUtil::createFolder(const std::string &folder_name){
    //std::cout<<"create folder "<<folder_name<<std::endl;
    boost::filesystem::path path(folder_name);
    boost::filesystem::create_directories(path);
  }

  void FileUtil::copyFile(const std::string &file_name,
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

  bool FileUtil::copyDir(boost::filesystem::path const & source, 
                         boost::filesystem::path const & destination){
    try{
      if(!fs::exists(source) || !fs::is_directory(source)){
        return false;
      }
      if(fs::exists(destination)){
        return false;
      }

      if(!fs::create_directory(destination)){
        std::cerr << "Unable to create destination directory"
            << destination.string() << std::endl;
        return false;
      }
    }
    catch(fs::filesystem_error const & e){
      std::cerr << e.what() << std::endl;
      return false;
    }

    for(fs::directory_iterator file(source); file != fs::directory_iterator(); ++file ){
      try{
        fs::path current(file->path());
        if(fs::is_directory(current)){
          if(!copyDir(current, destination / current.filename())){
            return false;
          }
        }
        else{
          fs::copy_file(current, destination / current.filename());
        }
      }
      catch(fs::filesystem_error const & e){
        std:: cerr << e.what() << std::endl;
      }
    }
    return true;
  }

  void FileUtil::delDir(const std::string &path){
    fs::path dir(path);

    if(fs::exists(dir))
      fs::remove_all(dir);
  }

  void FileUtil::delFile(const std::string &path){
    fs::path dir(path);

    if(fs::exists(dir))
      fs::remove(dir);
  }

  void FileUtil::cleanDir(const std::string &path){
    fs::path sp(path);
    fs::directory_iterator end_iter;
    std::string base = basename(absolute(sp).string());
    std::string output_table = "OUTPUT_TABLE", fasta = "fasta", msliagn = "msalign";

    for(fs::directory_iterator dir_iter(absolute(sp).parent_path()) ; dir_iter != end_iter ; ++dir_iter){
      std::string filename = dir_iter->path().string();
      if (filename.compare(0, base.length(), base) == 0){
        if (filename != absolute(sp).string()&&!fs::is_directory(dir_iter->path())&&
            filename.compare(filename.length() - output_table.length(), output_table.length(), output_table) != 0&&
            filename.compare(filename.length() - fasta.length(), fasta.length(), fasta) != 0&&
            filename.compare(filename.length() - msliagn.length(), msliagn.length(), msliagn) != 0){
          fs::remove(dir_iter->path());
        }
      } else if (filename.compare(filename.size() - 12, 12, "target_decoy") == 0){
        fs::remove(dir_iter->path());	
      } else if (filename.compare(filename.size() - 6, 6, "target") == 0){
        fs::remove(dir_iter->path());	
      }
    }
  }

  void clean_prefix(const fs::path & sp, const std::string & prefix) {
    fs::directory_iterator end_iter;
    for(fs::directory_iterator dir_iter(absolute(sp).parent_path()) ; dir_iter != end_iter ; ++dir_iter){
      std::string filename = dir_iter->path().string();
      if (filename.compare(0, prefix.length(), prefix) == 0){
        fs::remove(dir_iter->path());
      }
    }
  }

  void FileUtil::cleanDir(const std::string &fa_path, const std::string & sp_path) {
    fs::path fa(fa_path);
    fs::path sp(sp_path);
    std::string fa_base = absolute(fa).string();
    std::string sp_base = basename(absolute(sp).string());

    clean_prefix(fa, fa_base + "_");
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
    delFile(sp_base + ".CUTOFF_RESULT");            
    delFile(sp_base + ".LOCAL_RESULT");
    delFile(sp_base + ".EVALUE");
    delFile(sp_base + ".RAW_RESULT");      
    delFile(sp_base + ".FORMS");
    delFile(sp_base + ".FORM_RESULT");
    delFile(sp_base + ".FORM_FILTER_RESULT");

  }
}
