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

#include "topdia/common/topdia_para.hpp"

#include <iomanip>
#include <sstream>
#include <string>

#include "common/util/time_util.hpp"
#include "common/util/version.hpp"

namespace toppic {

namespace {

// Mirrors TopfdPara's printout so the two reports look identical: one label
// column (wide enough for the longest label) and fixed-width centered banners.
constexpr int kParaLabelWidth = 53;
constexpr int kParaBannerWidth = 55;

std::string banner(const std::string& prefix, const std::string& title) {
  int fill = kParaBannerWidth - 2 - static_cast<int>(title.size());
  if (fill < 2) fill = 2;
  int left = fill / 2;
  int right = fill - left;
  return prefix + std::string(left, '#') + " " + title + " " +
         std::string(right, '#');
}

}  // namespace

std::string TopdiaPara::getParaStr(const std::string& prefix,
                                   const std::string& sep,
                                   const TopfdParaPtr& topfd_para) const {
  std::stringstream output;
  const int w = kParaLabelWidth;
  auto kv = [&](const char* label) -> std::ostream& {
    return output << prefix << std::setw(w) << std::left << label << sep;
  };

  output << prefix << "TopDIA " << Version::getVersion() << std::endl;
  output << prefix << "Timestamp: " << time_util::getTimeStr() << std::endl;
  output << banner(prefix, "Parameters") << std::endl;

  output << topfd_para->getTopfdParaStr(prefix, sep);

  output << std::endl
         << banner(prefix, "TopDIA feature and pseudo-spectrum parameters")
         << std::endl;
  kv("MS2 Min scan number:") << topfd_para->getMs2MinScanNum() << std::endl;
  kv("MS1 ECScore cutoff:") << topfd_para->getMs1EcscoreCutoff() << std::endl;
  kv("MS2 ECScore cutoff:") << topfd_para->getMs2EcscoreCutoff() << std::endl;
  kv("Pseudo Score cutoff:") << pseudo_score_cutoff_ << std::endl;
  kv("Pseudo Min peak number:") << pseudo_min_peaks_ << std::endl;

  output << banner(prefix, "Parameters") << std::endl;
  return output.str();
}

}  // namespace toppic
