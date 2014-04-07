/*
 * folder_file.cpp
 *
 *  Created on: Apr 5, 2014
 *      Author: xunlikun
 */

#include <xpp/folder_file.hpp>

namespace prot {

FolderAndFile::FolderAndFile() {
  // TODO Auto-generated constructor stub

}

FolderAndFile::~FolderAndFile() {
  // TODO Auto-generated destructor stub
}


void createFolder(std::string folder_name){
  std::cout<<"create folder "<<folder_name<<std::endl;
  boost::filesystem::path path(folder_name);
  boost::filesystem::create_directories(path);
}

void copyFile(std::string file_name,std::string to_file,bool over_write){
  std::cout<<"copy file: "<<file_name<<" to "<<to_file<<std::endl;
  boost::filesystem::path from_path(file_name);
  boost::filesystem::path to_path(to_file);
  if(!boost::filesystem::exists(from_path)){
    std::cout<<"source file "<<file_name<<" not exist!"<<std::endl;
    return;
  }
  if(boost::filesystem::exists(to_path)){
    std::cout<<"target file "<<to_file<<" exist!"<<std::endl;
    if(over_write){
      std::cout<<"over writer the file :"<<to_file<<std::endl;
      boost::filesystem::remove(to_path);
    }
    else{
      return;
    }
  }


  boost::filesystem::copy(from_path,to_path);
}

} /* namespace prot */
