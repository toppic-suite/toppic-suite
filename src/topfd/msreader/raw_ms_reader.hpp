//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef TOPPIC_TOPFD_MSREADER_RAW_MS_READER_HPP_
#define TOPPIC_TOPFD_MSREADER_RAW_MS_READER_HPP_

#include "ms/spec/deconv_ms.hpp"
#include "ms/mzml/mzml_ms.hpp"
#include "ms/mzml/pw_ms_reader.hpp"

namespace toppic {

class RawMsReader {
 public:
  RawMsReader(const std::string & file_name, double isolation_window);

  RawMsReader(const std::string & file_name, 
              const std::string & activation, 
              double isolation_window);

  MzmlMsPtr getNextMs(double max_mass, int max_charge);

  int getInputSpNum() {return reader_ptr_->getInputSpNum();}

  void getMs1Peaks(PeakPtrVec2D &raw_peaks, double voltage);

  void getMs1Map(DeconvMsPtrVec &ms1_ptr_vec, DeconvMsPtrVec &ms2_ptr_vec, 
                 PeakPtrVec2D &ms1_raw_peaks, std::vector<double> &ms2_prec_mzs); 

  //void refinePrecChrg(MzmlMsPtr ms_one, MzmlMsPtr ms_two, 
  //                    double max_mass, int max_charge);

 private:
  PwMsReaderPtr reader_ptr_;
  MzmlMsPtr ms_one_; 

  //bool do_refine_prec_mass_ = true;

};

typedef std::shared_ptr<RawMsReader> RawMsReaderPtr;

}

#endif
