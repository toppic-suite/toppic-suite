/*
 * folder_file.hpp
 *
 *  Created on: Apr 5, 2014
 *      Author: xunlikun
 */

#ifndef FOLDER_FILE_HPP_
#define FOLDER_FILE_HPP_

#include <string>
#include <iostream>
#include <boost/filesystem.hpp>

namespace prot {

class FolderAndFile {
 public:
  FolderAndFile();
  virtual ~FolderAndFile();
};

void createFolder(std::string folder_name);
void copyFile(std::string file_name,std::string path,bool over_write);
} /* namespace prot */

#endif /* FOLDER_FILE_HPP_ */
