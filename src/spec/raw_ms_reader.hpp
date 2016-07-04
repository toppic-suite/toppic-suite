#ifndef PROT_SPEC_RAW_MS_READER_HPP_
#define PROT_SPEC_RAW_MS_READER_HPP_

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
  void readNext();

 private:
  std::string file_name_;
  int input_sp_num_;
  int input_sp_id_;
  int output_sp_id_;
  PeakPtrVec peak_list_;
  MsHeaderPtr header_ptr_;
  
  //pwiz reader
  pwiz::msdata::DefaultReaderList readers_;
	MSDataFilePtr msd_ptr_; 
  pwiz::msdata::SpectrumListPtr spec_list_ptr_;
   
};


}

#endif
