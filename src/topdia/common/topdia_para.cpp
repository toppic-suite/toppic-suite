//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include <iostream>
#include <sstream>
#include <iomanip>
#include "common/util/version.hpp"
#include "common/util/time_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/util/version.hpp"

#include "topdia/common/topdia_para.hpp"

namespace toppic {

std::string TopdiaPara::getParaStr(const std::string &prefix,
                                   const std::string &sep,
                                   TopfdParaPtr topfd_para) {
    std::stringstream output;
    int gap = 25;
    output << prefix << "TopDIA " << Version::getVersion() << std::endl;
    output << prefix << "Timestamp: " << time_util::getTimeStr() << std::endl;
    output << prefix << "###################### Parameters ######################" << std::endl;
    output << topfd_para->getTopfdParaStr(prefix, sep, gap); 
    output << prefix << std::setw(gap) << std::left
           << "MS2 Min scan number:        " << sep << topfd_para->getMs2MinScanNum() << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "MS1 ECScore cutoff:         " << sep  << topfd_para->getMs1EcscoreCutoff() << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "MS2 ECScore cutoff:         " << sep  << topfd_para->getMs2EcscoreCutoff() << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "Pseudo Score cutoff:        " << sep  << pseudo_score_cutoff_ << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "Pseudo Min peak number:     " << sep << pseudo_min_peaks_ << std::endl;
    output << prefix << std::setw(gap) << std::left
           << "Version:                    " << sep << Version::getVersion() << std::endl;
    output << prefix << "###################### Parameters ######################" << std::endl;
    return output.str();
  }

}  // namespace toppic
