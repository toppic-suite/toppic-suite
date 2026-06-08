//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include "topdia/common/topdia_para.hpp"

#include <iomanip>
#include <sstream>
#include <string>

#include "common/util/time_util.hpp"
#include "common/util/version.hpp"

namespace toppic {

std::string TopdiaPara::getParaStr(const std::string &prefix,
                                   const std::string &sep,
                                   const TopfdParaPtr &topfd_para) const {
    std::stringstream output;
    // Same label-column width as TopfdPara's printout so the TopDIA-specific
    // lines below line up with the TopFD parameters printed by getTopfdParaStr.
    const int gap = 53;
    output << prefix << "TopDIA " << Version::getVersion() << std::endl;
    output << prefix << "Timestamp: " << time_util::getTimeStr() << std::endl;
    output << prefix << "###################### Parameters ######################" << std::endl;
    output << topfd_para->getTopfdParaStr(prefix, sep);
    output << prefix << std::setw(gap) << std::left
           << "MS2 Min scan number:" << sep << topfd_para->getMs2MinScanNum() << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "MS1 ECScore cutoff:" << sep << topfd_para->getMs1EcscoreCutoff() << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "MS2 ECScore cutoff:" << sep << topfd_para->getMs2EcscoreCutoff() << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "Pseudo Score cutoff:" << sep << pseudo_score_cutoff_ << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "Pseudo Min peak number:" << sep << pseudo_min_peaks_ << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "Version:" << sep << Version::getVersion() << std::endl;
    output << prefix << "###################### Parameters ######################" << std::endl;
    return output.str();
  }

}  // namespace toppic
