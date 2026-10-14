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

#ifndef TOPPIC_FILTER_TAG_TAG_FILTER_MERGE_HPP_
#define TOPPIC_FILTER_TAG_TAG_FILTER_MERGE_HPP_

#include <string>

#include "filter/tag/tag_filter_mng.hpp"

namespace toppic {

namespace tag_filter_merge {

// Joins the tag filter's candidates to a mass filter's: for the var and
// one stages one file per proteoform type, <base>.<filter_ext>_<TYPE> and
// <base>.<tag_ext>_<TYPE> into <base>.<out_ext>_<TYPE>; for the other
// stages the three files without a type suffix. The filter's candidates of
// a spectrum come first, unchanged, then the tag candidates that are new:
//   one    a protein listed already with the same precursor variant gets
//          the tag shifts added to its N- and C-terminal shift lists;
//   var    a candidate of the same protein, range, N-terminal shift and
//          truncation N-terminal modification is not repeated;
//   multi, topmg   a protein listed with the same precursor variant is not
//          repeated.
void process(const std::string& sp_file_name, const std::string& filter_ext,
             const std::string& tag_ext, const std::string& out_ext,
             TagFilterMng::Stage stage);

}  // namespace tag_filter_merge

}  // namespace toppic

#endif
