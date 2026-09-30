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

#ifndef TOPPIC_TOPLIB_LIB_SEARCH_HPP_
#define TOPPIC_TOPLIB_LIB_SEARCH_HPP_

#include <map>
#include <string>
#include <vector>

#include "toplib/lib_spectrum.hpp"

namespace toppic {

namespace toplib {

// the dtype pandas' read_sql_query gives a column: int64 (integers only),
// float64 (numbers, or integers with a NULL) or object (text, or NULL only)
enum class PandasType { kInt, kFloat, kObject };

struct LibSearchColumn {
  PandasType type = PandasType::kObject;
  // one value per library row
  std::vector<LibValue> values;
};

// The rows of a library table (target_decoy_spectra_representatives) with a
// proteoform, targets and decoys, in table order.
struct LibSearchLibrary {
  int row_num = 0;
  // rows in the table, identified or not
  int table_row_num = 0;
  std::vector<double> precursor_masses;
  std::vector<double> precursor_charges;
  std::vector<double> flags;
  // kLibPeakNum values per row (row-major), NaN where a row has fewer peaks
  std::vector<double> masses;
  std::vector<double> norm_intensities;
  std::vector<double> charges;
  // the columns the result table copies, by name
  std::map<std::string, LibSearchColumn> columns;
};

enum class LibTolType { kPpm, kDa };

struct LibSearchParams {
  LibTolType tol_type = LibTolType::kPpm;
  // ppm of the library precursor mass, or Da
  double precursor_tol = 10;
  double fragment_ppm = 10;
  bool charge_match = false;
};

// a query spectrum's best library spectrum
struct LibMatch {
  // index into the query spectra (LibData::spectra)
  int query = 0;
  // index into the library rows
  int lib_row = 0;
  double distance = 0;
};

struct LibSearchResult {
  // the query spectra with a precursor intensity and at least two masses
  int searched_num = 0;
  // by query spectrum id
  std::vector<LibMatch> matches;
  // the distance cutoff whose FDR is closest to 1%
  bool cutoff_found = false;
  int cutoff_index = -1;
  double cutoff = 0;
  int cutoff_target_num = 0;
  int cutoff_decoy_num = 0;
  double cutoff_fdr = 0;
  // the matches written: the targets at the cutoff with a similarity of at
  // least 0.3, indices into matches
  std::vector<int> reported;
};

namespace lib_search {

// the library table read as ms_library_query.py's get_lib_representatives
// reads it
LibSearchLibrary readLibrary(const std::string& db_file_name);

// ms_query and fdr_cal
LibSearchResult search(const LibSearchLibrary& lib, const LibData& queries,
                       const LibSearchParams& params);

// query_res.tsv as pandas' to_csv writes it
void writeResult(const std::string& file_name, const LibSearchLibrary& lib,
                 const LibData& queries, const LibSearchResult& result);

// every match before the FDR step (spectrum_id, file_name, lib_index,
// representative_id, flag, calculated_cosine_dist), by query spectrum id
void writeMatches(const std::string& file_name, const LibSearchLibrary& lib,
                  const LibData& queries, const LibSearchResult& result);

}  // namespace lib_search

}  // namespace toplib

}  // namespace toppic

#endif
