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

#ifndef TOPPIC_TOPFD_MSREADER_RAW_MS_GROUP_FAIME_READER_HPP_
#define TOPPIC_TOPFD_MSREADER_RAW_MS_GROUP_FAIME_READER_HPP_

#include "ms/spec/raw_ms_group.hpp"
#include "topfd/msreader/pw_ms_reader.hpp"
#include "ms/env/match_env.hpp"

namespace toppic {

class RawMsGroupFaimeReader {
 public:
  RawMsGroupFaimeReader(const std::string & file_name, bool missing_level_one, 
                        std::string activation, double isolation_window,
                        int fraction_id);

  RawMsPtr readNextRawMs();

  RawMsGroupPtr getNextMsGroupPtrWithFaime();

  static void obtainPrecEnvs(RawMsGroupPtr ms_group_ptr, 
                             MatchEnvPtrVec &env_ptr_vec,
                             double max_mass,
                             int max_charge);

  int getInputSpNum() {return reader_ptr_->getInputSpNum();}
  PwMsReaderPtr getReaderPtr() {return reader_ptr_;}
  
 private:
  PwMsReaderPtr reader_ptr_;

  std::map<int, RawMsPtr> ms_one_ptr_map_;
  std::map<int, int> ms_one_scan_idx_map_;
  std::map<int, RawMsPtrVec> ms_two_ptr_vec_map_;
  std::map<int, int> last_ms_two_scan_map_;
  std::vector<int> ms_one_scans_;

  int total_ms_one_num_;
  int cur_ms_one_idx_;
  std::map<int, int> cur_last_ms_two_scan_map_;

  int fraction_id_;

  bool do_refine_prec_mass_ = true;
  bool missing_level_one_ = false;

  std::string activation_;
};

typedef std::shared_ptr<RawMsGroupFaimeReader> RawMsGroupFaimeReaderPtr;

}

#endif
