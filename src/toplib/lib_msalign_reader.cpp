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

#include "toplib/lib_msalign_reader.hpp"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "toplib/lib_numeric.hpp"

namespace toppic {

namespace toplib {

namespace lib_msalign_reader {

namespace {

// ------------------------------------------------------------ parsing
double toDouble(const std::string& s, const std::string& what) {
  char* end = nullptr;
  double v = std::strtod(s.c_str(), &end);
  if (end == s.c_str()) {
    LOG_ERROR("Cannot read the number \"" << s << "\" of " << what << ".");
    exit(EXIT_FAILURE);
  }
  return v;
}

long long toLong(const std::string& s, const std::string& what) {
  char* end = nullptr;
  long long v = std::strtoll(s.c_str(), &end, 10);
  if (end == s.c_str()) {
    LOG_ERROR("Cannot read the integer \"" << s << "\" of " << what << ".");
    exit(EXIT_FAILURE);
  }
  return v;
}

// the part of a field before its first ':' (the first of several precursors)
std::string firstValue(const std::string& v) { return v.substr(0, v.find(':')); }

// the text between the first and the second '=' (Python's line.split('=')[1])
std::string valueOf(const std::string& line) {
  size_t p = line.find('=');
  if (p == std::string::npos) {
    return "";
  }
  size_t q = line.find('=', p + 1);
  return line.substr(p + 1, q == std::string::npos ? std::string::npos : q - p - 1);
}

// the tab-separated columns of a peak line (whitespace if it has no tab)
std::vector<std::string> peakColumns(const std::string& line) {
  std::vector<std::string> cols = str_util::split(line, "\t");
  if (cols.size() < 3) {
    cols.clear();
    std::string cur;
    for (char c : line) {
      if (std::isspace(static_cast<unsigned char>(c))) {
        if (!cur.empty()) {
          cols.push_back(cur);
          cur.clear();
        }
      } else {
        cur += c;
      }
    }
    if (!cur.empty()) {
      cols.push_back(cur);
    }
  }
  return cols;
}

// pandas' read_csv NA spellings
bool isNa(const std::string& s) {
  static const std::set<std::string> na = {
      "",     "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN",
      "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL",
      "NaN",  "None", "n/a", "nan", "null"};
  return na.count(s) > 0;
}

// a table cell pandas reads as a number
bool isNumber(const std::string& s) {
  double v;
  return lib_numeric::pandasToDouble(s, v);
}

// The identification columns in the order of TopLib's library table, with
// the PrSM table columns each can come from: TopLib's name first (TopPIC
// 1.7's table), then the names TopPIC 1.9 and 2.0 gave the same column.
struct IdColumnSource {
  std::string name;
  std::vector<std::string> headers;
};

const std::vector<IdColumnSource>& idColumnMap() {
  static const std::vector<IdColumnSource> m = {
      {"proteoform_id", {"Proteoform ID"}},
      {"protein_accession", {"Protein accession"}},
      {"protein_description", {"Protein description"}},
      {"first_residue", {"First residue"}},
      {"last_residue", {"Last residue"}},
      {"protein_sequence", {"Database protein sequence"}},
      {"proteoform", {"Proteoform"}},
      {"e_value", {"E-value"}},
      {"spectrum_level_q_value", {"Spectrum-level Q-value"}},
      {"proteoform_level_q_value", {"Proteoform-level Q-value"}},
      {"num_unexpected_modifications",
       {"#unexpected modifications", "# Unexpected modifications"}},
      {"unexpected_modifications",
       {"unexpected modifications", "Unexpected modifications"}},
      {"fixed_ptms", {"Fixed PTMs"}},
      {"adjusted_precursor_mass", {"Adjusted precursor mass"}},
      {"prsm_id", {"Prsm ID"}},
      {"fragmentation", {"Fragmentation"}},
      {"num_original_peaks", {"#peaks", "# Masses"}},
      {"proteoform_intensity", {"Proteoform intensity"}},
      {"feature_intensity", {"Feature intensity"}},
      {"feature_score", {"Feature score"}},
      {"feature_apex_time", {"Feature apex time"}},
      {"num_protein_hits", {"#Protein hits", "# Protein hits"}},
      {"proteoform_mass", {"Proteoform mass"}},
      {"protein_n_terminal_form", {"Protein N-terminal form"}},
      {"num_variable_ptms", {"#variable PTMs", "# Variable PTMs"}},
      {"variable_ptms", {"variable PTMs", "Variable PTMs"}},
      {"num_matched_peaks",
       {"#matched peaks", "# Matched masses", "# Matched exp fragments"}},
      {"num_matched_fragment_ions",
       {"#matched fragment ions", "# Matched fragments",
        "# Matched theo fragments"}},
      {"special_amion_acids", {"Special amino acids"}},
      {"miscore", {"MIScore"}}};
  return m;
}

// the columns ms_library_building.py reads unconditionally
const int kRequiredIdColumnNum = 10;

// Turns text columns into typed ones as pandas does: a column whose non-NA
// values in all rows of the table (table_rows, matched to a spectrum or
// not) are numbers is REAL (pandas' float64 after the .loc assignment into
// the spectrum frame, integers included; the values parsed as pandas parses
// them, which is not always correctly rounded), otherwise TEXT; an all-NA
// column is REAL. raw[s][c] is the text of column c for spectrum s, empty rows for
// unidentified spectra.
void typeIdColumns(LibData& data,
                   const std::vector<std::vector<std::string>>& table_rows,
                   const std::vector<std::vector<std::string>>& raw,
                   const std::vector<bool>& identified) {
  size_t col_num = data.id_columns.size();
  for (size_t c = 0; c < col_num; c++) {
    bool numeric = true;
    for (const std::vector<std::string>& row : table_rows) {
      if (!isNa(row[c]) && !isNumber(row[c])) {
        numeric = false;
        break;
      }
    }
    data.id_columns[c].type = numeric ? LibColType::kReal : LibColType::kText;
  }
  for (size_t s = 0; s < raw.size(); s++) {
    if (!identified[s]) {
      continue;
    }
    std::vector<LibValue>& values = data.spectra[s].id_values;
    values.assign(col_num, LibValue());
    for (size_t c = 0; c < col_num; c++) {
      const std::string& v = raw[s][c];
      if (isNa(v)) {
        continue;
      }
      if (data.id_columns[c].type == LibColType::kReal) {
        double d = 0;
        lib_numeric::pandasToDouble(v, d);
        values[c] = d;
      } else {
        values[c] = v;
      }
    }
  }
}

// the NaN bookkeeping of the optional msalign fields
void noteOptional(LibData& data, const LibSpectrum& sp) {
  data.any_ms_level |= sp.ms_level.has_value();
  data.all_ms_level &= sp.ms_level.has_value();
  data.any_ms_one_id |= sp.ms_one_id.has_value();
  data.all_ms_one_id &= sp.ms_one_id.has_value();
  data.any_ms_one_scan |= sp.ms_one_scan.has_value();
  data.all_ms_one_scan &= sp.ms_one_scan.has_value();
  data.any_window_begin |= sp.precursor_window_begin.has_value();
  data.any_window_end |= sp.precursor_window_end.has_value();
  data.any_precursor_mz |= sp.precursor_mz.has_value();
}

// the blocks between BEGIN IONS and END IONS lines, each a list of lines
std::vector<std::vector<std::string>> readBlocks(const std::string& file_name,
                                                 bool substring) {
  std::ifstream in(file_name);
  if (!in) {
    LOG_ERROR("Cannot open the msalign file " << file_name << ".");
    exit(EXIT_FAILURE);
  }
  std::vector<std::vector<std::string>> blocks;
  std::vector<std::string> cur;
  bool inside = false;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    bool begin = substring ? line.find("BEGIN IONS") != std::string::npos
                           : line == "BEGIN IONS";
    bool end = substring ? line.find("END IONS") != std::string::npos
                         : line == "END IONS";
    if (begin) {
      cur.clear();
      inside = true;
    }
    if (inside) {
      cur.push_back(line);
    }
    if (end && inside) {
      blocks.push_back(cur);
      inside = false;
    }
  }
  return blocks;
}

bool isMassLine(const std::string& line) {
  return !line.empty() && std::isdigit(static_cast<unsigned char>(line[0]));
}

// the masses of a block no larger than the precursor mass, in file order
void readPeaks(const std::vector<std::string>& block, LibSpectrum& sp) {
  for (const std::string& line : block) {
    if (!isMassLine(line)) {
      continue;
    }
    std::vector<std::string> cols = peakColumns(line);
    if (cols.size() < 3) {
      LOG_ERROR("Cannot read the peak line \"" << line << "\".");
      exit(EXIT_FAILURE);
    }
    double mass = toDouble(cols[0], "a fragment mass");
    double inte = toDouble(cols[1], "a fragment intensity");
    int charge = static_cast<int>(toLong(cols[2], "a fragment charge"));
    if (mass > sp.precursor_mass) {
      continue;
    }
    sp.masses.push_back(mass);
    sp.intensities.push_back(inte);
    sp.charges.push_back(charge);
  }
}

// the precursor intensities of a colon-separated list
std::vector<double> intensityList(const std::string& v) {
  std::vector<double> result;
  for (const std::string& s : str_util::split(v, ":")) {
    result.push_back(toDouble(s, "a precursor intensity"));
  }
  return result;
}

}  // namespace

// ------------------------------------------------------------ TopFD + table
LibData readMsalignAndPrsms(const std::string& msalign_file_name,
                            const std::string& prsm_file_name,
                            bool python_compatible) {
  LibData data;
  // read_msalign: every field is the first header line containing its key
  for (const std::vector<std::string>& block :
       readBlocks(msalign_file_name, true)) {
    bool has_mass = false;
    for (const std::string& line : block) {
      has_mass |= isMassLine(line);
    }
    if (!has_mass) {
      continue;
    }
    auto find = [&block](const std::string& key) -> const std::string* {
      for (const std::string& line : block) {
        if (line.find(key) != std::string::npos) {
          return &line;
        }
      }
      return nullptr;
    };
    auto need = [&find, &msalign_file_name](const std::string& key) {
      const std::string* line = find(key);
      if (line == nullptr) {
        LOG_ERROR("A spectrum in " << msalign_file_name << " has no " << key
                                   << " field.");
        exit(EXIT_FAILURE);
      }
      return valueOf(*line);
    };
    LibSpectrum sp;
    sp.file_name = need("FILE_NAME");
    sp.precursor_mass = toDouble(firstValue(need("PRECURSOR_MASS")),
                                 "PRECURSOR_MASS");
    std::string inte_value = need("PRECURSOR_INTENSITY");
    sp.precursor_intensity =
        toDouble(firstValue(inte_value), "PRECURSOR_INTENSITY");
    sp.precursor_charge = static_cast<int>(
        toLong(firstValue(need("PRECURSOR_CHARGE")), "PRECURSOR_CHARGE"));
    sp.precursor_feature_id =
        toLong(firstValue(need("PRECURSOR_FEATURE_ID")), "PRECURSOR_FEATURE_ID");
    sp.spectrum_id = toLong(need("SPECTRUM_ID="), "SPECTRUM_ID");
    sp.scan = toLong(need("SCANS="), "SCANS");
    sp.retention_time = toDouble(need("RETENTION_TIME"), "RETENTION_TIME");
    sp.precursor_intensities = intensityList(inte_value);
    if (const std::string* l = find("PRECURSOR_MZ")) {
      sp.precursor_mz = toDouble(firstValue(valueOf(*l)), "PRECURSOR_MZ");
    }
    if (const std::string* l = find("TITLE")) {
      sp.title = valueOf(*l);
    }
    if (const std::string* l = find("LEVEL")) {
      sp.ms_level = toLong(valueOf(*l), "LEVEL");
    }
    if (const std::string* l = find("MS_ONE_ID")) {
      sp.ms_one_id = toLong(valueOf(*l), "MS_ONE_ID");
    }
    if (const std::string* l = find("MS_ONE_SCAN")) {
      sp.ms_one_scan = toLong(valueOf(*l), "MS_ONE_SCAN");
    }
    if (const std::string* l = find("PRECURSOR_WINDOW_BEGIN")) {
      sp.precursor_window_begin =
          toDouble(valueOf(*l), "PRECURSOR_WINDOW_BEGIN");
    }
    if (const std::string* l = find("PRECURSOR_WINDOW_END")) {
      sp.precursor_window_end = toDouble(valueOf(*l), "PRECURSOR_WINDOW_END");
    }
    if (const std::string* l = find("ACTIVATION")) {
      sp.activation = valueOf(*l);
    }
    // skip a spectrum with a zero precursor mass
    if (sp.precursor_mass == 0.0) {
      continue;
    }
    readPeaks(block, sp);
    noteOptional(data, sp);
    data.spectra.push_back(std::move(sp));
  }

  // the PrSM table (tab-separated): the header is the first line with a
  // Spectrum ID column; the lines before it (TopPIC's parameter block,
  // which TopLib removes beforehand) are skipped
  std::ifstream in(prsm_file_name);
  if (!in) {
    LOG_ERROR("Cannot open the PrSM file " << prsm_file_name << ".");
    exit(EXIT_FAILURE);
  }
  std::string line;
  std::vector<std::string> header;
  std::map<std::string, int> pos;
  while (pos.count("Spectrum ID") == 0 && std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    header = str_util::split(line, "\t");
    pos.clear();
    for (size_t i = 0; i < header.size(); i++) {
      pos[header[i]] = static_cast<int>(i);
    }
  }
  if (pos.count("Spectrum ID") == 0) {
    LOG_ERROR("The PrSM file " << prsm_file_name << " has no Spectrum ID column.");
    exit(EXIT_FAILURE);
  }
  std::vector<int> src_col;
  const auto& cmap = idColumnMap();
  for (size_t i = 0; i < cmap.size(); i++) {
    // the first of the column's names in the table (TopLib's name only in
    // the Python-compatible mode)
    size_t name_num = python_compatible ? 1 : cmap[i].headers.size();
    int col = -1;
    for (size_t k = 0; k < name_num && col < 0; k++) {
      auto it = pos.find(cmap[i].headers[k]);
      if (it != pos.end()) {
        col = it->second;
      }
    }
    if (col < 0) {
      if (static_cast<int>(i) < kRequiredIdColumnNum) {
        LOG_ERROR("The PrSM file " << prsm_file_name << " has no "
                                   << cmap[i].headers[0] << " column.");
        exit(EXIT_FAILURE);
      }
      continue;
    }
    data.id_columns.push_back({cmap[i].name, LibColType::kReal});
    src_col.push_back(col);
  }
  std::multimap<long long, size_t> by_id;
  for (size_t s = 0; s < data.spectra.size(); s++) {
    by_id.insert({data.spectra[s].spectrum_id, s});
  }
  std::vector<std::vector<std::string>> raw(data.spectra.size());
  std::vector<bool> identified(data.spectra.size(), false);
  std::vector<std::vector<std::string>> table_rows;
  int spec_col = pos["Spectrum ID"];
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty()) {
      continue;
    }
    std::vector<std::string> cols = str_util::split(line, "\t");
    cols.resize(header.size());
    long long spec_id = toLong(cols[spec_col], "Spectrum ID");
    std::vector<std::string> values;
    for (int c : src_col) {
      values.push_back(cols[c]);
    }
    auto range = by_id.equal_range(spec_id);
    for (auto it = range.first; it != range.second; ++it) {
      // a later row of the same spectrum replaces an earlier one
      raw[it->second] = values;
      identified[it->second] = true;
    }
    table_rows.push_back(std::move(values));
  }
  typeIdColumns(data, table_rows, raw, identified);
  return data;
}

// ------------------------------------------------------------ library msalign
LibData readLibraryMsalign(const std::string& msalign_file_name) {
  LibData data;
  const std::vector<std::pair<std::string, std::string>> id_keys = {
      {"proteoform_id", ""},
      {"protein_accession", "PROTEIN_ACCESSION"},
      {"protein_description", ""},
      {"first_residue", "FIRST_RESIDUE_POSITION"},
      {"last_residue", ""},
      {"protein_sequence", "DATABASE_SEQUENCE"},
      {"proteoform", "PROTEOFORM"},
      {"e_value", "E_VALUE"},
      {"spectrum_level_q_value", ""},
      {"proteoform_level_q_value", ""},
      {"unexpected_modifications", "UNEXPECTED_MODIFICATIONS"},
      {"fixed_ptms", "FIXED_MODIFICATIONS"}};
  for (const auto& k : id_keys) {
    data.id_columns.push_back({k.first, LibColType::kReal});
  }
  std::map<std::pair<std::string, std::string>, long long> form_ids;
  std::vector<std::vector<std::string>> raw;
  std::vector<bool> identified;
  // the identifications of all blocks (the converted table's rows)
  std::vector<std::vector<std::string>> table_rows;
  for (const std::vector<std::string>& block :
       readBlocks(msalign_file_name, false)) {
    std::map<std::string, std::string> head;
    bool has_mass = false;
    for (const std::string& line : block) {
      if (isMassLine(line)) {
        has_mass = true;
      } else {
        size_t p = line.find('=');
        if (p != std::string::npos) {
          head[line.substr(0, p)] = line.substr(p + 1);
        }
      }
    }
    auto get = [&head](const std::string& key) -> const std::string* {
      auto it = head.find(key);
      return it == head.end() ? nullptr : &it->second;
    };
    auto need = [&get, &msalign_file_name](const std::string& key) {
      const std::string* v = get(key);
      if (v == nullptr) {
        LOG_ERROR("A spectrum in " << msalign_file_name << " has no " << key
                                   << " field.");
        exit(EXIT_FAILURE);
      }
      return *v;
    };
    // the proteoform ids count every identified block in file order
    std::vector<std::string> values;
    bool ided = get("PROTEOFORM") != nullptr && !get("PROTEOFORM")->empty();
    if (ided) {
      std::string acc = get("PROTEIN_ACCESSION") ? *get("PROTEIN_ACCESSION") : "";
      const std::string& form = *get("PROTEOFORM");
      auto key = std::make_pair(acc, form);
      if (form_ids.count(key) == 0) {
        long long next = static_cast<long long>(form_ids.size());
        form_ids[key] = next;
      }
      // residues: the capital letters outside [...]
      long long residues = 0;
      int depth = 0;
      for (char c : form) {
        if (c == '[') {
          depth++;
        } else if (c == ']') {
          depth--;
        } else if (depth == 0 && c >= 'A' && c <= 'Z') {
          residues++;
        }
      }
      long long first = toLong(need("FIRST_RESIDUE_POSITION"),
                               "FIRST_RESIDUE_POSITION");
      for (const auto& k : id_keys) {
        if (k.first == "proteoform_id") {
          values.push_back(std::to_string(form_ids[key]));
        } else if (k.first == "last_residue") {
          values.push_back(std::to_string(first + residues - 1));
        } else if (k.second.empty()) {
          values.push_back("");
        } else {
          const std::string* v = get(k.second);
          values.push_back(v == nullptr ? "" : *v);
        }
      }
      table_rows.push_back(values);
    }
    if (!has_mass) {
      continue;
    }
    LibSpectrum sp;
    sp.file_name = need("MZML_FILE_NAME");
    sp.precursor_mass = toDouble(firstValue(need("PRECURSOR_MONOISOTOPIC_MASS")),
                                 "PRECURSOR_MONOISOTOPIC_MASS");
    std::string inte_value = need("PRECURSOR_INTENSITY");
    sp.precursor_intensity =
        toDouble(firstValue(inte_value), "PRECURSOR_INTENSITY");
    sp.precursor_charge = static_cast<int>(
        toLong(firstValue(need("PRECURSOR_CHARGE")), "PRECURSOR_CHARGE"));
    sp.precursor_feature_id =
        toLong(firstValue(need("PRECURSOR_FEATURE_ID")), "PRECURSOR_FEATURE_ID");
    sp.spectrum_id = toLong(need("MS2_ID"), "MS2_ID");
    sp.scan = toLong(need("MS2_SCAN"), "MS2_SCAN");
    sp.retention_time = toDouble(need("MS2_RETENTION_TIME"), "MS2_RETENTION_TIME");
    sp.precursor_intensities = intensityList(inte_value);
    if (const std::string* v = get("PRECURSOR_MONOISOTOPIC_MZ")) {
      sp.precursor_mz = toDouble(firstValue(*v), "PRECURSOR_MONOISOTOPIC_MZ");
    }
    if (const std::string* v = get("TITLE")) {
      sp.title = *v;
    }
    if (const std::string* v = get("LEVEL")) {
      sp.ms_level = toLong(*v, "LEVEL");
    }
    if (const std::string* v = get("MS1_ID")) {
      sp.ms_one_id = toLong(*v, "MS1_ID");
    }
    if (const std::string* v = get("MS1_SCAN")) {
      sp.ms_one_scan = toLong(*v, "MS1_SCAN");
    }
    if (const std::string* v = get("PRECURSOR_WINDOW_BEGIN")) {
      sp.precursor_window_begin = toDouble(*v, "PRECURSOR_WINDOW_BEGIN");
    }
    if (const std::string* v = get("PRECURSOR_WINDOW_END")) {
      sp.precursor_window_end = toDouble(*v, "PRECURSOR_WINDOW_END");
    }
    if (const std::string* v = get("ACTIVATION")) {
      sp.activation = *v;
    }
    if (sp.precursor_mass == 0.0) {
      continue;
    }
    readPeaks(block, sp);
    noteOptional(data, sp);
    data.spectra.push_back(std::move(sp));
    raw.push_back(values);
    identified.push_back(ided);
  }
  typeIdColumns(data, table_rows, raw, identified);
  return data;
}

}  // namespace lib_msalign_reader

}  // namespace toplib

}  // namespace toppic
