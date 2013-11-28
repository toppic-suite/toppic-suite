#ifndef PROT_MSALIGN_READER_HPP_
#define PROT_MSALIGN_READER_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include "deconv_ms.hpp"

namespace prot {

class MsAlignReader {
 public:
  MsAlignReader(const char *spectrum_file, ActivationPtrVec activation_list);

  std::vector<std::string> readOneSpectrum();

  void readNext();

  DeconvMsPtr getNextMs();

  std::vector<std::string> getSpectrumStr() {
    return spectrum_str_;
  }

  void close();

 private:
  std::ifstream input_;
  ActivationPtrVec activation_list_;
  std::vector<std::string> spectrum_str_;
  int current_ = 1;
  DeconvMsPtr deconv_ms_ptr_ = DeconvMsPtr(nullptr);
};

int countSpNum(const char *spectrum_file);

}
#endif
