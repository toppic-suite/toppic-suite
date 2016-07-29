#ifndef PROT_FEATURE_RAW_MS_READER_HPP_
#define PROT_FEATURE_RAW_MS_READER_HPP_

#include <memory>
#include <vector>
#include <string>

#include "pwiz/data/msdata/DefaultReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"

#include "spec/raw_ms.hpp"

namespace prot {

typedef std::shared_ptr<pwiz::msdata::MSDataFile> MSDataFilePtr;

class RawMsReader {
 public:
  RawMsReader(std::string &file_name);
  int readNext();
  PeakPtrVec getPeakList() {return peak_list_;}
  MsHeaderPtr getHeaderPtr() {return header_ptr_;}
  int getInputSpNum() {return input_sp_num_;}

 private:
  std::string file_name_;
  int input_sp_num_;
  int input_sp_id_;
  int output_sp_id_;

  int ms1_cnt = 0;
  int ms2_cnt = 0;
  PeakPtrVec peak_list_;
  MsHeaderPtr header_ptr_;
  
  //pwiz reader
  pwiz::msdata::DefaultReaderList readers_;
	MSDataFilePtr msd_ptr_; 
  pwiz::msdata::SpectrumListPtr spec_list_ptr_;
};

typedef std::shared_ptr<RawMsReader> RawMsReaderPtr;

}

#endif
