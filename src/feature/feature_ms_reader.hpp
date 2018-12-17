//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_MS_READER_HPP_
#define PROT_FEATURE_MS_READER_HPP_

#include <memory>
#include <vector>

#include "spec/raw_ms.hpp"
#include "feature/raw_ms_reader.hpp"

namespace toppic {

class FeatureMsReader {
 public:
  FeatureMsReader(const std::string & file_name);

  RawMsPtr getNextMs(double prec_win_size, int max_charge);

  void refinePrecChrg(RawMsPtr ms_one, RawMsPtr ms_two, 
                      double prec_win_size, int max_charge);

  int getInputSpNum() {return reader_ptr_->getInputSpNum();}

 private:
  RawMsReaderPtr reader_ptr_;
  RawMsPtr ms_one_; 

  bool do_refine_prec_mass_ = true;

};

typedef std::shared_ptr<FeatureMsReader> FeatureMsReaderPtr;

}

#endif
