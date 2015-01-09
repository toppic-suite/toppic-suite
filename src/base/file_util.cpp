#include <iostream>
#include <sstream>

#include "base/logger.hpp"
#include "base/file_util.hpp"

#ifdef PROT_LINUX

#include <unistd.h>

#endif

namespace fs = boost::filesystem;

namespace prot
{

#ifdef PROT_LINUX
//linux
std::string getExecutiveDir(const std::string &argv_0){
	
  int buffer_size = 1024;
  char* buffer = new char[buffer_size];
  size_t len = readlink ("/proc/self/exe", buffer, buffer_size);
  std::string file_name;
  if (len >= 0){
    buffer[len] = '\0';
    file_name = std::string (buffer);
  }
  else{
    LOG_ERROR("Can not find executive directory!");
    return file_name = "";
  }
  delete[] buffer;
  fs::path full_path(file_name);
  std::string exe_dir = full_path.remove_filename().string();
  return exe_dir;
}

#else

//windows
std::string getExecutiveDir(const std::string &argv_0){
	
  fs::path full_path;
  full_path = fs::system_complete( fs::path(argv_0));

  std::string exe_dir = full_path.remove_filename().string();
  return exe_dir;
}

#endif

std::string basename(const std::string &s){
  size_t dot_pos = s.find_last_of(".");
  if (dot_pos < s.length()){
    return s.substr(0, dot_pos);
  }
  return s;
}

std::string directory(const std::string &s){
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

bool copyDir(boost::filesystem::path const & source, boost::filesystem::path const & destination){
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

void delDir(const std::string &path){
  fs::path dir(path);
	
  if(fs::exists(dir))
    fs::remove_all(dir);
}

void delFile(const std::string &path){
  fs::path dir(path);
	
  if(fs::exists(dir))
    fs::remove(dir);
}

void cleanDir(const std::string &path){
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

}
