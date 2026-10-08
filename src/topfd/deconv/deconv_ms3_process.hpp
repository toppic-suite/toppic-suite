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

#ifndef TOPPIC_TOPFD_DECONV_DECONV_MS3_PROCESS_HPP_
#define TOPPIC_TOPFD_DECONV_DECONV_MS3_PROCESS_HPP_

#include <map>
#include <memory>
#include <string>

#include "ms/spec/deconv_ms.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

// The deconvolution of MS3 spectra (the dissociation of an MS2 fragment
// ion; toppic_test experiment 146), a stage after the MS/MS deconvolution
// that reads the mzML file again for the MS3 scans without a proton
// transfer stage (those are PTCR spectra, handled before). The precursor of
// an MS3 spectrum is the species of its parent MS2 scan's deconvolution
// (<base>_ms2.msalign) whose isotopic envelope has a peak in the MS3
// isolation window, the most intense when several; without one, the
// selected ion of the file at every charge from 1 to the MS2 precursor's
// (the file's charge state is not usable). The spectrum is deconvoluted with
// the MS/MS method bounded by the precursor's mass and charge and written
// to <base>_ms3.msalign (LEVEL=3, MS_ONE_* the cycle's MS1 scan, MS_TWO_*
// the parent MS2 scan and its precursor); <base>_ms3.tsv lists one row per
// MS3 scan.
class DeconvMs3Process {
 public:
  explicit DeconvMs3Process(const TopfdParaPtr& topfd_para_ptr);

  void process();

 private:
  // the MS1 / MS2 spectra of the msalign file: scan number -> spectrum
  void readMsAlign(const std::string& file_name,
                   std::map<int, DeconvMsPtr>& ms_map);

  TopfdParaPtr topfd_para_ptr_;
};

using DeconvMs3ProcessPtr = std::shared_ptr<DeconvMs3Process>;

}  // namespace toppic

#endif
