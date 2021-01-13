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

#ifndef TOPPIC_TOPFD_MSREADER_RAW_MS_GROUP_READER_HPP_
#define TOPPIC_TOPFD_MSREADER_RAW_MS_GROUP_READER_HPP_

#include "ms/spec/raw_ms_group.hpp"
#include "topfd/msreader/pw_ms_reader.hpp"
#include "ms/env/match_env.hpp"

namespace toppic {

class RawMsGroupReader {
 public:
  RawMsGroupReader(const std::string & file_name, bool missing_level_one, std::string activation, 
                   int fraction_id);

  RawMsPtr readNextRawMs();

  RawMsGroupPtr getNextMsGroupPtr();

  static void obtainPrecEnvs(RawMsGroupPtr ms_group_ptr, 
                             MatchEnvPtrVec &env_ptr_vec,
                             double prec_win_size, int max_charge);

  int getInputSpNum() {return reader_ptr_->getInputSpNum();}

 private:
  PwMsReaderPtr reader_ptr_;
  RawMsPtr ms_one_ptr_; 
  int fraction_id_;

  bool do_refine_prec_mass_ = true;
  bool missing_level_one_ = false;

  std::string activation_;
};

typedef std::shared_ptr<RawMsGroupReader> RawMsGroupReaderPtr;

}

#endif
