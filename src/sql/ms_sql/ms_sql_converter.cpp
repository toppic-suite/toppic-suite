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

#include "sql/ms_sql/ms_sql_converter.hpp"

#include <cstddef>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "pwiz/data/msdata/DefaultReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"

#include "common/util/file_util.hpp"
#include "sql/ms_sql/ms_sql_writer.hpp"

namespace toppic {

namespace ms_sql_converter {

std::string getDbFileName(const std::string& spec_file_name) {
  return file_util::basename(spec_file_name) + "_3d.db";
}

void convert(const std::string& spec_file_name, double mz_size,
             double rt_divider) {
  std::cout << "Reading the MS1 peaks of " << spec_file_name << " - started."
            << std::endl;
  pwiz::msdata::DefaultReaderList readers;
  pwiz::msdata::MSDataFile msd(spec_file_name, &readers);
  pwiz::msdata::SpectrumListPtr spec_list_ptr = msd.run.spectrumListPtr;

  std::vector<MsSqlPeak> peaks;
  int ms1_scan_num = 0;
  size_t spec_num = spec_list_ptr->size();
  for (size_t i = 0; i < spec_num; i++) {
    // Read the header first; the peak list is only decoded for MS1 scans.
    pwiz::msdata::SpectrumPtr spec_ptr = spec_list_ptr->spectrum(i, false);
    pwiz::msdata::SpectrumInfo spec_info(*spec_ptr);
    if (spec_info.msLevel != 1) {
      continue;
    }
    ms1_scan_num++;
    int rt = static_cast<int>(spec_info.retentionTime * 1000);  // ms
    spec_ptr = spec_list_ptr->spectrum(i, true);
    std::vector<pwiz::msdata::MZIntensityPair> pairs;
    spec_ptr->getMZIntensityPairs(pairs);
    for (const pwiz::msdata::MZIntensityPair& pair : pairs) {
      if (pair.intensity <= 0.0) {
        continue;
      }
      MsSqlPeak peak;
      peak.mz = pair.mz;
      peak.inte = pair.intensity;
      peak.rt = rt;
      peaks.push_back(peak);
    }
  }
  std::cout << "Reading the MS1 peaks of " << spec_file_name << " - finished: "
            << peaks.size() << " peaks in " << ms1_scan_num << " MS1 scans."
            << std::endl;

  std::string db_file_name = getDbFileName(spec_file_name);
  std::cout << "Writing " << db_file_name << " - started." << std::endl;
  MsSqlWriter writer(db_file_name);
  writer.write(std::move(peaks), ms1_scan_num, mz_size, rt_divider);
  std::cout << "Writing " << db_file_name << " - finished." << std::endl;
}

}  // namespace ms_sql_converter

}  // namespace toppic
