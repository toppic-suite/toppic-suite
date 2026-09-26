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

#ifndef TOPPIC_FILTER_MNG_PTM_DIAG_INDEX_HPP_
#define TOPPIC_FILTER_MNG_PTM_DIAG_INDEX_HPP_

namespace toppic {

// The PRM diagonal index file that two filters share: the one-shift filter
// (TopPIC's one-shift search, TopMG's ASF-One PTM filter) and the
// multiple-shift filter (TopPIC's two-shift search, TopMG's ASF-Diagonal
// filter) build the same MassMatch from a block's proteoforms, so topindex
// writes it once, as the diagonal file of the one-shift family, and both
// filters read it. The two managers therefore take the file name and the
// column limits it is built at from here, so that they cannot diverge; the
// filters also check a file's column number against their limit when they
// read it.
namespace ptm_diag_index {

constexpr const char* kFileName = "ptm_diag_index";

// the column limits: the largest proteoform mass (Da) and the scale of the
// integer masses
constexpr int kMaxProteoformMass = 50000;
constexpr int kFilterScale = 100;

}  // namespace ptm_diag_index

}  // namespace toppic

#endif
