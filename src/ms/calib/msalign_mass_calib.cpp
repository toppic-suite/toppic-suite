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

#include "ms/calib/msalign_mass_calib.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "common/util/logger.hpp"

namespace toppic {

namespace msalign_mass_calib {

namespace {

std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  size_t start = 0;
  while (true) {
    size_t pos = s.find(sep, start);
    if (pos == std::string::npos) {
      out.push_back(s.substr(start));
      return out;
    }
    out.push_back(s.substr(start, pos - start));
    start = pos + 1;
  }
}

// the values of a ':'-separated header list, empty fields left out
std::vector<std::string> headerValues(const std::string& line) {
  std::vector<std::string> out;
  size_t eq = line.find('=');
  if (eq == std::string::npos) return out;
  for (const std::string& v : split(line.substr(eq + 1), ':')) {
    if (!v.empty()) out.push_back(v);
  }
  return out;
}

std::string fmt6(double v) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.6f", v);
  return buf;
}

std::string join(const std::vector<std::string>& v) {
  std::string s;
  for (size_t i = 0; i < v.size(); i++) {
    if (i > 0) s += ":";
    s += v[i];
  }
  return s;
}

bool startsWith(const std::string& s, const char* prefix) {
  return s.rfind(prefix, 0) == 0;
}

// the corrected lines of one spectrum (BEGIN IONS ... END IONS)
std::vector<std::string> correctBlock(const MassCalibModel& model,
                                      const std::vector<std::string>& block) {
  double proton = mass_constant::getProtonMass();
  bool has_rt = false;
  double rt = 0.0;
  std::vector<std::string> mzs, charges, masses;
  for (const std::string& line : block) {
    if (startsWith(line, "RETENTION_TIME=")) {
      rt = std::atof(line.c_str() + 15);
      has_rt = true;
    } else if (startsWith(line, "PRECURSOR_MZ=")) {
      mzs = headerValues(line);
    } else if (startsWith(line, "PRECURSOR_CHARGE=")) {
      charges = headerValues(line);
    } else if (startsWith(line, "PRECURSOR_MASS=")) {
      masses = headerValues(line);
    }
  }
  if (!has_rt) {
    LOG_ERROR("A spectrum of the msalign file has no RETENTION_TIME!");
    std::exit(EXIT_FAILURE);
  }
  std::vector<std::string> new_mzs, new_masses;
  size_t prec_num = std::min(mzs.size(), std::min(charges.size(), masses.size()));
  for (size_t i = 0; i < prec_num; i++) {
    double mz = std::atof(mzs[i].c_str());
    int z = std::atoi(charges[i].c_str());
    double m = std::atof(masses[i].c_str());
    double e = model.getErrorPpm(1, mz, rt) * 1e-6;
    new_mzs.push_back(fmt6(mz * (1.0 - e)));
    new_masses.push_back(fmt6(m - e * (m + z * proton)));
  }
  std::vector<std::string> out;
  out.reserve(block.size());
  for (const std::string& line : block) {
    if (startsWith(line, "PRECURSOR_MZ=") && !new_mzs.empty()) {
      out.push_back("PRECURSOR_MZ=" + join(new_mzs));
    } else if (startsWith(line, "PRECURSOR_MASS=") && !new_masses.empty()) {
      out.push_back("PRECURSOR_MASS=" + join(new_masses));
    } else if (!line.empty() &&
               std::isdigit(static_cast<unsigned char>(line[0]))) {
      std::vector<std::string> f = split(line, '\t');
      if (f.size() >= 3) {
        double m = std::atof(f[0].c_str());
        int z = std::atoi(f[2].c_str());
        if (z > 0) {
          double e = model.getErrorPpm(2, (m + z * proton) / z, rt) * 1e-6;
          f[0] = fmt6(m - e * (m + z * proton));
        }
        std::string s = f[0];
        for (size_t i = 1; i < f.size(); i++) s += "\t" + f[i];
        out.push_back(s);
      } else {
        out.push_back(line);
      }
    } else {
      out.push_back(line);
    }
  }
  return out;
}

}  // namespace

int process(const MassCalibModel& model, const std::string& in_file_name,
            const std::string& out_file_name) {
  std::ifstream in(in_file_name);
  if (!in) {
    LOG_ERROR("Cannot open the msalign file " << in_file_name << "!");
    std::exit(EXIT_FAILURE);
  }
  std::ofstream out(out_file_name);
  if (!out) {
    LOG_ERROR("Cannot write the msalign file " << out_file_name << "!");
    std::exit(EXIT_FAILURE);
  }
  int sp_num = 0;
  bool in_block = false;
  std::vector<std::string> block;
  std::string line;
  while (std::getline(in, line)) {
    if (startsWith(line, "BEGIN IONS")) {
      block.clear();
      block.push_back(line);
      in_block = true;
    } else if (in_block) {
      block.push_back(line);
      if (startsWith(line, "END IONS")) {
        for (const std::string& l : correctBlock(model, block)) out << l << "\n";
        in_block = false;
        sp_num++;
      }
    } else {
      out << line << "\n";
    }
  }
  // an unterminated last spectrum is copied unchanged
  if (in_block) {
    for (const std::string& l : block) out << l << "\n";
  }
  return sp_num;
}

}  // namespace msalign_mass_calib

}  // namespace toppic
