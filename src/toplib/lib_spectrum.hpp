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

#ifndef TOPPIC_TOPLIB_LIB_SPECTRUM_HPP_
#define TOPPIC_TOPLIB_LIB_SPECTRUM_HPP_

#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace toppic {

namespace toplib {

// A cell of the library table: NULL (TopLib's pandas NaN), an integer, a
// real or a text value.
using LibValue = std::variant<std::monostate, long long, double, std::string>;

// The SQLite column types TopLib's pandas data frames are written with:
// int64 -> INTEGER, float64 -> REAL, object -> TEXT.
enum class LibColType { kInteger, kReal, kText };

// A column of identification data (from the PrSM table or from the headers
// of a library msalign file), named as in the library table.
struct LibIdColumn {
  std::string name;
  LibColType type = LibColType::kReal;
};

// One MS/MS spectrum as TopLib's read_msalign reads it, with the
// identification get_spectra_by_file attaches to it.
struct LibSpectrum {
  std::string file_name;
  long long spectrum_id = 0;
  std::string title;
  long long scan = 0;
  double retention_time = 0;
  std::optional<long long> ms_level;
  std::optional<long long> ms_one_id;
  std::optional<long long> ms_one_scan;
  std::optional<double> precursor_window_begin;
  std::optional<double> precursor_window_end;
  std::string activation;
  std::optional<double> precursor_mz;
  // the first precursor of the spectrum
  double precursor_mass = 0;
  double precursor_intensity = 0;
  int precursor_charge = 0;
  long long precursor_feature_id = 0;
  // the intensities of all its precursors (the PPIR filter)
  std::vector<double> precursor_intensities;
  // the fragment masses no larger than the precursor mass, in file order
  std::vector<double> masses;
  std::vector<double> intensities;
  std::vector<int> charges;
  // one value per LibData::id_columns; empty when unidentified
  std::vector<LibValue> id_values;
};

// The spectra of one input and the columns that describe them.
struct LibData {
  std::vector<LibSpectrum> spectra;
  std::vector<LibIdColumn> id_columns;
  // TopLib drops a column that is NaN for every spectrum and stores an
  // integer column with a NaN as REAL: whether any / every spectrum has
  // the optional msalign fields
  bool any_ms_level = false, all_ms_level = true;
  bool any_ms_one_id = false, all_ms_one_id = true;
  bool any_ms_one_scan = false, all_ms_one_scan = true;
  bool any_window_begin = false, any_window_end = false;
  bool any_precursor_mz = false;
};

}  // namespace toplib

}  // namespace toppic

#endif
