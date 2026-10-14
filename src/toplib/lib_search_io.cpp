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

#include "toplib/lib_search.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <system_error>
#include <variant>
#include <vector>

#include <sqlite3.h>

#include "common/util/logger.hpp"
#include "toplib/lib_builder.hpp"
#include "toplib/lib_numeric.hpp"

namespace toppic {

namespace toplib {

namespace lib_search {

namespace {

const char kTableName[] = "target_decoy_spectra_representatives";

// the columns get_lib_representatives selects besides the peaks; the
// script fails without any of them
const std::vector<std::string>& selectedColumns() {
  static const std::vector<std::string> names = {
      "representative_id", "file_id", "file_name", "spectrum_id", "scan",
      "precursor_mass", "precursor_intensity", "precursor_charge",
      "retention_time", "precursor_feature_id", "proteoform_id",
      "protein_accession", "protein_description", "first_residue",
      "last_residue", "protein_sequence", "proteoform", "e_value",
      "spectrum_level_q_value", "proteoform_level_q_value", "flag"};
  return names;
}

// the library columns of query_res.tsv, each written as "<name>(lib)" after
// representative_id and similarity
const std::vector<std::string>& resultLibColumns() {
  static const std::vector<std::string> names = {
      "file_name", "spectrum_id", "scan", "retention_time", "proteoform_id",
      "protein_accession", "precursor_mass", "precursor_charge",
      "precursor_feature_id", "protein_description", "protein_sequence",
      "first_residue", "last_residue", "proteoform", "e_value",
      "spectrum_level_q_value", "proteoform_level_q_value"};
  return names;
}

// a value as Python's sqlite3 returns it
LibValue columnValue(sqlite3_stmt* stmt, int c) {
  switch (sqlite3_column_type(stmt, c)) {
    case SQLITE_INTEGER:
      return static_cast<long long>(sqlite3_column_int64(stmt, c));
    case SQLITE_FLOAT:
      return sqlite3_column_double(stmt, c);
    case SQLITE_NULL:
      return LibValue();
    default: {
      const unsigned char* text = sqlite3_column_text(stmt, c);
      int n = sqlite3_column_bytes(stmt, c);
      return std::string(reinterpret_cast<const char*>(text), n);
    }
  }
}

// DataFrame.from_records(coerce_float=True), over all rows of the table
PandasType pandasType(const std::vector<LibValue>& values) {
  bool any_int = false;
  bool any_float = false;
  bool any_null = false;
  bool any_text = false;
  for (const LibValue& v : values) {
    any_null |= std::holds_alternative<std::monostate>(v);
    any_int |= std::holds_alternative<long long>(v);
    any_float |= std::holds_alternative<double>(v);
    any_text |= std::holds_alternative<std::string>(v);
  }
  if (any_text || (!any_int && !any_float)) {
    return PandasType::kObject;
  }
  if (any_float || any_null) {
    return PandasType::kFloat;
  }
  return PandasType::kInt;
}

// a value of a float64 column (NaN for NULL)
double numberOf(const LibValue& v, const std::string& column) {
  if (const long long* i = std::get_if<long long>(&v)) {
    return static_cast<double>(*i);
  }
  if (const double* d = std::get_if<double>(&v)) {
    return *d;
  }
  if (std::holds_alternative<std::string>(v)) {
    LOG_ERROR("The library column " << column << " has the text \""
                                    << std::get<std::string>(v)
                                    << "\" where a number is expected.");
    exit(EXIT_FAILURE);
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// csv.QUOTE_MINIMAL with a tab delimiter
std::string csvField(const std::string& s) {
  if (s.find_first_of("\t\"\n\r") == std::string::npos) {
    return s;
  }
  std::string q = "\"";
  for (char c : s) {
    if (c == '"') {
      q += "\"\"";
    } else {
      q += c;
    }
  }
  return q + "\"";
}

std::string floatCell(double v) {
  return std::isnan(v) ? "" : lib_numeric::pythonRepr(v);
}

// a cell as pandas' to_csv writes a value of a column of this dtype: int64
// as an integer, float64 as Python's repr (an empty cell for NaN), object
// as str() of the Python value (an empty cell for None)
std::string cellText(const LibValue& v, PandasType type,
                     const std::string& column) {
  if (std::holds_alternative<std::monostate>(v)) {
    return "";
  }
  if (type == PandasType::kFloat) {
    return floatCell(numberOf(v, column));
  }
  if (const long long* i = std::get_if<long long>(&v)) {
    return std::to_string(*i);
  }
  if (const double* d = std::get_if<double>(&v)) {
    return lib_numeric::pythonRepr(*d);
  }
  return csvField(std::get<std::string>(v));
}

std::string libCell(const LibSearchLibrary& lib, const std::string& column,
                    int row) {
  const LibSearchColumn& col = lib.columns.at(column);
  return cellText(col.values[row], col.type, column);
}

std::ofstream openOutput(const std::string& file_name) {
  std::filesystem::path parent =
      std::filesystem::path(file_name).parent_path();
  if (!parent.empty()) {
    std::error_code ec;
    std::filesystem::create_directories(parent, ec);
  }
  std::ofstream out(file_name);
  if (!out) {
    LOG_ERROR("Cannot write " << file_name << ".");
    exit(EXIT_FAILURE);
  }
  return out;
}

}  // namespace

LibSearchLibrary readLibrary(const std::string& db_file_name) {
  sqlite3* db = nullptr;
  if (sqlite3_open_v2(db_file_name.c_str(), &db, SQLITE_OPEN_READONLY,
                      nullptr) != SQLITE_OK) {
    LOG_ERROR("Cannot open the library " << db_file_name << ".");
    exit(EXIT_FAILURE);
  }
  std::string sql = std::string("SELECT * FROM \"") + kTableName + "\"";
  sqlite3_stmt* stmt = nullptr;
  if (sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr) != SQLITE_OK) {
    LOG_ERROR("The library " << db_file_name << " has no table " << kTableName
                             << ": " << sqlite3_errmsg(db));
    exit(EXIT_FAILURE);
  }
  int col_num = sqlite3_column_count(stmt);
  std::map<std::string, int> col_index;
  for (int c = 0; c < col_num; c++) {
    col_index.emplace(sqlite3_column_name(stmt, c), c);
  }
  std::vector<std::vector<LibValue>> table(col_num);
  int rc;
  while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
    for (int c = 0; c < col_num; c++) {
      table[c].push_back(columnValue(stmt, c));
    }
  }
  if (rc != SQLITE_DONE) {
    LOG_ERROR("Cannot read the library " << db_file_name << ": "
                                         << sqlite3_errmsg(db));
    exit(EXIT_FAILURE);
  }
  sqlite3_finalize(stmt);
  sqlite3_close(db);

  auto need = [&col_index, &db_file_name](const std::string& name) {
    auto it = col_index.find(name);
    if (it == col_index.end()) {
      LOG_ERROR("The library " << db_file_name << " has no column " << name
                               << ".");
      exit(EXIT_FAILURE);
    }
    return it->second;
  };
  for (const std::string& name : selectedColumns()) {
    need(name);
  }
  std::vector<int> mass_cols, inte_cols, charge_cols, norm_cols;
  for (int k = 0; k < kLibPeakNum; k++) {
    std::string n = std::to_string(k);
    mass_cols.push_back(need("mass_" + n));
    inte_cols.push_back(need("intensity_" + n));
    charge_cols.push_back(need("charge_" + n));
    norm_cols.push_back(need("norm_intensity_" + n));
  }

  // dropna(subset='proteoform'): the identified rows, targets and decoys
  LibSearchLibrary lib;
  const std::vector<LibValue>& proteoforms = table[need("proteoform")];
  lib.table_row_num = static_cast<int>(proteoforms.size());
  std::vector<int> rows;
  for (int r = 0; r < lib.table_row_num; r++) {
    if (!std::holds_alternative<std::monostate>(proteoforms[r])) {
      rows.push_back(r);
    }
  }
  lib.row_num = static_cast<int>(rows.size());

  auto numbers = [&table, &rows, &need](const std::string& name) {
    const std::vector<LibValue>& col = table[need(name)];
    std::vector<double> v;
    v.reserve(rows.size());
    for (int r : rows) {
      v.push_back(numberOf(col[r], name));
    }
    return v;
  };
  lib.precursor_masses = numbers("precursor_mass");
  lib.precursor_charges = numbers("precursor_charge");
  lib.flags = numbers("flag");

  auto peaks = [&table, &rows](const std::vector<int>& cols) {
    std::vector<double> v(rows.size() * kLibPeakNum);
    for (size_t i = 0; i < rows.size(); i++) {
      for (int k = 0; k < kLibPeakNum; k++) {
        v[i * kLibPeakNum + k] =
            numberOf(table[cols[k]][rows[i]], "of the peaks");
      }
    }
    return v;
  };
  lib.masses = peaks(mass_cols);
  lib.norm_intensities = peaks(norm_cols);
  lib.charges = peaks(charge_cols);

  // the columns copied into the output keep the dtype of the whole table
  std::vector<std::string> copied = resultLibColumns();
  copied.push_back("representative_id");
  copied.push_back("flag");
  for (const std::string& name : copied) {
    const std::vector<LibValue>& values = table[need(name)];
    LibSearchColumn col;
    col.type = pandasType(values);
    for (int r : rows) {
      col.values.push_back(values[r]);
    }
    lib.columns[name] = std::move(col);
  }
  return lib;
}

void writeResult(const std::string& file_name, const LibSearchLibrary& lib,
                 const LibData& queries, const LibSearchResult& result) {
  std::ofstream out = openOutput(file_name);
  out << "spectrum_id\tfile_name\trepresentative_id\tsimilarity";
  for (const std::string& name : resultLibColumns()) {
    out << "\t" << name << "(lib)";
  }
  out << "\n";
  for (int m : result.reported) {
    const LibMatch& match = result.matches[m];
    const LibSpectrum& sp = queries.spectra[match.query];
    out << sp.spectrum_id << "\t" << csvField(sp.file_name) << "\t"
        << libCell(lib, "representative_id", match.lib_row) << "\t"
        << floatCell(1 - match.distance);
    for (const std::string& name : resultLibColumns()) {
      out << "\t" << libCell(lib, name, match.lib_row);
    }
    out << "\n";
  }
}

void writeMatches(const std::string& file_name, const LibSearchLibrary& lib,
                  const LibData& queries, const LibSearchResult& result) {
  std::ofstream out = openOutput(file_name);
  out << "spectrum_id\tfile_name\tlib_index\trepresentative_id\tflag\t"
      << "calculated_cosine_dist\n";
  for (const LibMatch& match : result.matches) {
    const LibSpectrum& sp = queries.spectra[match.query];
    out << sp.spectrum_id << "\t" << csvField(sp.file_name) << "\t"
        << match.lib_row << "\t"
        << libCell(lib, "representative_id", match.lib_row) << "\t"
        << libCell(lib, "flag", match.lib_row) << "\t"
        << floatCell(match.distance) << "\n";
  }
}

}  // namespace lib_search

}  // namespace toplib

}  // namespace toppic
