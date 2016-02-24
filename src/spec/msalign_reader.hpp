#ifndef PROT_SPEC_MSALIGN_READER_HPP_
#define PROT_SPEC_MSALIGN_READER_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"

namespace prot {

class MsAlignReader {
 public:
  MsAlignReader(const std::string &file_name, int group_spec_num, 
                ActivationPtr activation_ptr);

  std::vector<std::string> readOneSpectrum();

  void readNext();

  DeconvMsPtr getNextMs();

  SpectrumSetPtr getNextSpectrumSet(SpParaPtr sp_para_ptr);

  void close();

 private:
  std::string file_name_;
  int group_spec_num_;
  std::ifstream input_;
  std::vector<std::string> spectrum_str_vec_;
  int current_ = 0;
  ActivationPtr activation_ptr_;
  DeconvMsPtr deconv_ms_ptr_ = DeconvMsPtr(nullptr);

};

typedef std::shared_ptr<MsAlignReader> MsAlignReaderPtr;

}
#endif
