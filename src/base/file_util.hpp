#ifndef PROT_BASE_FILE_UTIL_HPP_
#define PROT_BASE_FILE_UTIL_HPP_

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

#include <string>

namespace prot {

class FileUtil {
 public:
  static std::string getFileSeparator();

  static std::string getExecutiveDir(const std::string &argv_0);

  static std::string basename(const std::string &s);

  static std::string directory(const std::string &s);

  static void createFolder(const std::string &folder_name);

  static void copyFile(const std::string &file_name, const std::string &path,
                       bool over_write);

  static bool copyDir(boost::filesystem::path const & source, 
                      boost::filesystem::path const & destination);

  static void delDir(const std::string &path);

  static void delFile(const std::string &path);

  static void cleanDir(const std::string &path);
};

}

#endif
