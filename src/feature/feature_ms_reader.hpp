#ifndef PROT_FEATURE_MS_READER_HPP_
#define PROT_FEATURE_MS_READER_HPP_

#include <memory>
#include <vector>

#include "spec/raw_ms.hpp"
#include "spec/raw_ms_reader.hpp"

namespace prot {

class FeatureMsReader {
 public:
  FeatureMsReader(std::string &file_name);
  
  RawMsPtr getNextMs(double prec_win_size);

  void refinePrecChrg(RawMsPtr ms_one, RawMsPtr ms_two, 
                      double prec_win_size);

  int getInputSpNum() {return reader_ptr_->getInputSpNum();}

 private:
  RawMsReaderPtr reader_ptr_;
  RawMsPtr ms_one_; 

  bool do_refine_prec_mass_ = true;

};

typedef std::shared_ptr<FeatureMsReader> FeatureMsReaderPtr;

}

#endif
