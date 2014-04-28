#include <iostream>

#include <sstream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "base/file_util.hpp"

namespace fs = boost::filesystem;

namespace prot {

std::string getExeDir(const std::string &argv_0) {
  fs::path full_path;
  full_path = fs::system_complete( fs::path(argv_0));

  std::string exe_dir = full_path.remove_filename().string();
  return exe_dir;
}

}
