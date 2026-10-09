// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "filter/tag/tag_filter_mng.hpp"

#include <algorithm>

#include "common/util/str_util.hpp"

namespace toppic {

TagFilterMng::TagFilterMng(const PrsmParaPtr& prsm_para_ptr, Stage stage,
                           int thread_num, const std::string& output_file_ext)
    : prsm_para_ptr_(prsm_para_ptr),
      stage_(stage),
      thread_num_(thread_num),
      output_file_ext_(output_file_ext) {}

const std::vector<std::string>& TagFilterMng::getStageNames() {
  static const std::vector<std::string> names = {"var", "one", "multi",
                                                 "topmg"};
  return names;
}

const std::string& TagFilterMng::getStageName(Stage stage) {
  return getStageNames()[static_cast<int>(stage)];
}

bool TagFilterMng::parseStages(const std::string& value,
                               std::set<Stage>* stages) {
  stages->clear();
  if (value == "") {
    return true;
  }
  const std::vector<std::string>& names = getStageNames();
  std::vector<std::string> parts = str_util::split(value, ",");
  for (size_t i = 0; i < parts.size(); i++) {
    std::string part = parts[i];
    str_util::trim(part);
    if (part == "all") {
      for (size_t k = 0; k < names.size(); k++) {
        stages->insert(static_cast<Stage>(k));
      }
      continue;
    }
    size_t pos = std::find(names.begin(), names.end(), part) - names.begin();
    if (pos == names.size()) {
      stages->clear();
      return false;
    }
    stages->insert(static_cast<Stage>(pos));
  }
  return true;
}

const std::vector<double>& TagFilterMng::getPrecErrorVec() const {
  SpParaPtr sp_para_ptr = prsm_para_ptr_->getSpParaPtr();
  switch (stage_) {
    case Stage::kVar:
      return sp_para_ptr->getVarPtmSearchPrecErrorVec();
    case Stage::kMulti:
      return sp_para_ptr->getMultiShiftSearchPrecErrorVec();
    case Stage::kOne:
    case Stage::kTopmg:
    default:
      // the hosted TopMG search filters with TopMG's one-shift filter
      return sp_para_ptr->getOneShiftSearchPrecErrorVec();
  }
}

bool TagFilterMng::searchSpectrum(int spec_id) const {
  if (skip_spec_ids_.count(spec_id) > 0) {
    return false;
  }
  return search_all_ || search_spec_ids_.count(spec_id) > 0;
}

}  // namespace toppic
