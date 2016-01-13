#include "base/logger.hpp"
#include "spec/msalign_util.hpp"

namespace prot {

int MsAlignUtil::countSpNum(const std::string &spectrum_file_name) {
  MsAlignReader reader(spectrum_file_name, 1);
  int cnt = 0;
  DeconvMsPtr deconv_ms_ptr;
  while ((deconv_ms_ptr = reader.getNextMs()) != nullptr) {
    cnt++;
  }
  reader.close();
  return cnt;
}

void MsAlignUtil::geneSpIndex(const std::string &spectrum_file_name) {
  int sp_num = countSpNum(spectrum_file_name); 
  std::ofstream index_output;
  std::string index_file_name = spectrum_file_name + "_index";
  index_output.open(index_file_name.c_str(), std::ios::out);
  index_output << sp_num << std::endl;
  index_output.close();
}

int MsAlignUtil::getSpNum(const std::string &spectrum_file_name) {
  std::ifstream index_input;
  std::string index_file_name = spectrum_file_name + "_index";
  index_input.open(index_file_name.c_str(), std::ios::in);
  std::string line;
  std::getline(index_input, line);
  int sp_num = std::stoi(line);
  LOG_DEBUG("Get sp number " << sp_num);
  return sp_num;
}

}
