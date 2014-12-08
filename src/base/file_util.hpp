#ifndef PROT_FILE_UTIL_HPP_
#define PROT_FILE_UTIL_HPP_

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

#include <string>

namespace prot
{

#define FILE_SEPARATOR "/"

std::string getExecutiveDir(const std::string &argv_0);

std::string basename(const std::string &s);

std::string directory(const std::string &s);

void createFolder(const std::string &folder_name);

void copyFile(const std::string &file_name, const std::string &path,
              bool over_write);

bool copyDir(boost::filesystem::path const & source, boost::filesystem::path const & destination);

void delDir(const std::string &path);

void delFile(const std::string &path);

void cleanDir(const std::string &path);

}

#endif
