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

#ifndef TOPPIC_PRSM_PRSM_POST_MASS_MATCH_HPP_
#define TOPPIC_PRSM_PRSM_POST_MASS_MATCH_HPP_

#include <string>

#include "para/prsm_para.hpp"

namespace toppic {

// Post mass matching (on by default in toppic; --disable-post-match turns it
// off): after the search, match the
// theoretical fragment masses of each PrSM that none of TopFD's deconvoluted
// masses matched against the centroided MS/MS peaks that TopFD stored in its
// SQLite database, using the method of MSPathFinderT (Informed-Proteomics):
// for every charge state from 1 to min(15, precursor charge + 2), the
// theoretical isotopic envelope of the fragment is placed in the spectrum and
// accepted if its most abundant isotopic peak is present, at least
// min_peak_num of its isotopic peaks are present, and the observed intensities
// have a Pearson correlation >= 0.7 or a Bhattacharyya distance <= 0.03 with
// the theoretical ones. Each accepted fragment adds one experimental
// monoisotopic mass to the spectrum, whose confidence score is the EnvCNN
// score of its envelope, computed as TopFD does for the deconvoluted masses.
namespace prsm_post_mass_match {

// Reads the PrSMs of <spectrum base name>.<input_file_ext> with their spectra,
// runs the post mass matching, and writes
//   <base>_post_ms2.msalign     the spectra with the added masses (the input
//                              msalign file is left untouched),
//   <base>_post_ms2.<output_file_ext>  the PrSMs with the recounted matched
//                              masses and fragments,
//   <base>.sqlite               the added masses as new rows of the ms2_env
//                              and ms2_env_peak tables (the spectrum file's
//                              database, sqlFileName),
//   <base>_post_ms2.feature    a copy of the TopFD feature file, so the
//                              downstream steps find it under the new name.
// Returns the name of the new spectrum file, which the downstream steps read
// instead of the input one. The SQLite database must exist (see
// sqlFileName); if it is missing the run stops with an error. The EnvCNN
// model must have been loaded with onnx_env_cnn::initModel.
// The database is the spectrum file's (sqlFileName): TopFD's, toppic's
// copy of it for -g, or the corrected copy for --mass-calib, whose peaks are
// read as they are.
std::string process(const PrsmParaPtr& prsm_para_ptr,
                    const std::string& input_file_ext,
                    const std::string& output_file_ext, int min_peak_num);

// The SQLite database that goes with an msalign spectrum file:
// <path>/<base>.sqlite, where <base> is the file name without its extension
// and a trailing "_ms2" (TopFD writes <base>_ms2.msalign and <base>.sqlite),
// and without "_resid" before it (the residual spectra of --second-precursor
// share the database of the spectra they are made from). toppic's copies of
// a spectrum file have their own copies of TopFD's database, so that a run
// on them never writes into TopFD's file: <base>_mono_ms2.msalign of -g
// (--average-mass) <base>_mono.sqlite, the corrected spectra of
// --mass-calib <base>_calib_ms2.msalign <base>_calib.sqlite and, after -g,
// <base>_mono_calib_ms2.msalign <base>_mono_calib.sqlite (written by
// TopPICProgress_multi_file and TopPIC_mass_calib). The argument parser
// checks that TopFD's database exists before the search.
std::string sqlFileName(const std::string& sp_file_name);

// TopFD's base name of a spectrum file, <path>/<base> for <base>_ms2.msalign
// and toppic's copies of it (the suffixes "_resid", "_calib" and "_mono"
// removed): the base of TopFD's database and MS1 spectra
std::string topfdBaseName(const std::string& sp_file_name);

// the base name of the post-matched files of a spectrum file: the
// database's base name with "_post_ms2", <base>_post_ms2,
// <base>_mono_post_ms2 for the copy of -g, <base>_calib_post_ms2 and
// <base>_mono_calib_post_ms2 for the corrected spectra of --mass-calib; the
// residual spectra of --second-precursor keep "_resid":
// <base>_resid_post_ms2, <base>_mono_calib_resid_post_ms2
std::string postBaseName(const std::string& sp_file_name);

}  // namespace prsm_post_mass_match

}  // namespace toppic

#endif
