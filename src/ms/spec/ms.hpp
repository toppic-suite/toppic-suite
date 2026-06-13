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

#ifndef TOPPIC_MS_SPEC_MS_HPP_
#define TOPPIC_MS_SPEC_MS_HPP_

#include <cstddef>
#include <utility>
#include <vector>

#include "ms/spec/ms_header.hpp"

namespace toppic {

template <class T>
class Ms {
 public:
  // A spectrum keeps its header and peak list, so both are sink parameters:
  // taken by value and moved into the members (see CLAUDE.md parameter
  // passing).
  Ms(MsHeaderPtr header_ptr, std::vector<T> peak_ptr_list)
      : header_ptr_(std::move(header_ptr)),
        peak_ptr_list_(std::move(peak_ptr_list)) {}

  MsHeaderPtr getMsHeaderPtr() const { return header_ptr_; }

  size_t size() const { return peak_ptr_list_.size(); }

  T getPeakPtr(int i) const { return peak_ptr_list_[i]; }

  const std::vector<T>& getPeakPtrVec() const { return peak_ptr_list_; }

  void setPeakPtrVec(std::vector<T> peak_ptr_list) {
    peak_ptr_list_ = std::move(peak_ptr_list);
  }

 private:
  MsHeaderPtr header_ptr_;
  std::vector<T> peak_ptr_list_;
};

}  // namespace toppic
#endif
