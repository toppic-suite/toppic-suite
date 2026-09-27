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

#ifndef TOPPIC_TOPFD_DECONV_FLASHDECONV_MS2_PROCESS_HPP_
#define TOPPIC_TOPFD_DECONV_FLASHDECONV_MS2_PROCESS_HPP_

#include <map>
#include <memory>
#include <string>

#include "ms/feature/spec_feature.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

// Deconvolutes the MS/MS spectra of an mzML file with the FLASHDeconv
// algorithm and writes the same output files as DeconvMs2Process: one msalign
// file per thread that is merged at the end, and, when the SQLite output is
// enabled, the spectra with their envelopes.
//
// The two classes differ only in how one spectrum is deconvoluted.
// DeconvMs2Process builds a candidate envelope for every peak and charge and
// picks a non-overlapping subset by dynamic programming (DeconvSingleSp);
// this class instead follows FLASHDeconv (OpenMS, src/topp/FLASHDeconv.cpp and
// ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm): the peaks are binned on a
// logarithmic m/z axis, masses supported by several charge states are
// harvested from the charge-mass matrix, and every candidate is scored with
// the isotope cosine score, the charge distribution score and QScore before
// the harmonic and charge-error artifacts are removed.
class FlashDeconvMs2Process {
 public:
  FlashDeconvMs2Process(const TopfdParaPtr& topfd_para_ptr,
                        const std::string& output_filename_ext);

  void process();

 private:
  TopfdParaPtr topfd_para_ptr_;
  std::string output_filename_ext_;

  void readSpecFeature(const std::string& feat_file_name,
                       std::map<int, SpecFeaturePtrVec>& feat_map);
};

using FlashDeconvMs2ProcessPtr = std::shared_ptr<FlashDeconvMs2Process>;

}  // namespace toppic

#endif
