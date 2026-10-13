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
//   <base>_post_ms2.feature    a copy of the TopFD feature file, so the
//                              downstream steps find it under the new name,
// and, when toppic_db_name is not empty, appends the added masses as new rows
// of the ms2_env and ms2_env_peak tables of toppic_db_name, toppic's database
// of the spectrum file (toppicSqlFileName of the file given to toppic;
// empty with --no-sql). Returns the name of the new spectrum file, which the
// downstream steps read instead of the input one. The peaks are read from
// TopFD's database topfd_db_name (sqlFileName of the file given to toppic),
// opened read only as immutable (no -wal / -shm files; it must not change
// during the run), in one pass over its ms2_peak table (toppic_test
// experiments 164 and 168). With --mass-calib, calib_model_file_name is the
// calibration model (<base>_mass_calib.json) and every peak's m/z is
// corrected with its MS/MS model, as in toppic's corrected copy of the
// database; otherwise it is empty. If TopFD's database is missing the run
// stops with an error. The EnvCNN model must have been loaded with
// onnx_env_cnn::initModel.
std::string process(const PrsmParaPtr& prsm_para_ptr,
                    const std::string& input_file_ext,
                    const std::string& output_file_ext, int min_peak_num,
                    const std::string& topfd_db_name,
                    const std::string& calib_model_file_name,
                    const std::string& toppic_db_name);

// TopFD's SQLite database of a spectrum file: <path>/<base>.sqlite, where
// <base> is the file name without its extension and a trailing "_ms2"
// (TopFD writes <base>_ms2.msalign and <base>.sqlite), and without "_resid"
// before it. toppic only reads it (the argument check, the copy of its
// envelope tables, or with --mass-calib of every table, into
// toppicSqlFileName, and the post mass matching); topmg writes its PrSM
// tables into it (not with --no-sql).
std::string sqlFileName(const std::string& sp_file_name);

// toppic's SQLite database of a spectrum file given to toppic,
// <path>/<base>_toppic.sqlite for <base>_ms2.msalign and TopFD's
// <base>.sqlite (toppic_test experiment 156): written again before each
// search of the file, with TopFD's tables ms2_env and ms2_env_peak (toppic_test
// experiment 168; every table before it) or, with --mass-calib when it fits
// a model, the corrected copy of every table of TopFD's database; the post
// mass matching appends its envelopes to it and the report writes toppic's
// tables into it. A reader takes a table it lacks from TopFD's database.
// -A, --mass-calib and --second-precursor use the database of the given
// file (TopPICProgress_multi_file passes it to every step).
std::string toppicSqlFileName(const std::string& sp_file_name);

// TopFD's base name of a spectrum file, <path>/<base> for <base>_ms2.msalign
// and toppic's copies of it (the suffixes "_resid", "_calib" and "_mono"
// removed): the base of TopFD's MS1 spectra
std::string topfdBaseName(const std::string& sp_file_name);

// the base name of the post-matched files of a spectrum file: the spectrum
// file's base name with "_post_ms2", <base>_post_ms2,
// <base>_mono_post_ms2 for the copy of -A, <base>_calib_post_ms2 and
// <base>_mono_calib_post_ms2 for the corrected spectra of --mass-calib; the
// residual spectra of --second-precursor keep "_resid":
// <base>_resid_post_ms2, <base>_mono_calib_resid_post_ms2
std::string postBaseName(const std::string& sp_file_name);

}  // namespace prsm_post_mass_match

}  // namespace toppic

#endif
