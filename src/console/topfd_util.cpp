//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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
#include "console/topfd_util.hpp"

namespace toppic {

namespace topfd_util {

std::string geneArgumentStr(std::map<std::string,std::string> &arguments,
                            const std::string & prefix) {
  std::stringstream output;
  output << prefix << "TopFD " << Version::getVersion() << std::endl;
  // TIME_STAMP_STR is replaced later
  output << prefix << "Timestamp: " << time_util::TIME_STAMP_STR << std::endl;
  output << prefix << "###################### Parameters ######################" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Data type: " << "Centroid" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Maximum charge: " << arguments["maxCharge"] << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Maximum monoisotopic mass: " << arguments["maxMass"] << " Dalton" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Error tolerance: " << arguments["mzError"] << " m/z" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "MS1 signal/noise ratio: " << arguments["msOneSnRatio"] << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "MS/MS signal/noise ratio: " << arguments["msTwoSnRatio"] << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Precursor window size: " << arguments["precWindow"] << " m/z" << std::endl;
  //output << prefix << std::setw(40) << std::left 
  //    << "Do final filtering: " << para_ptr->do_final_filtering_ << std::endl;
  output << prefix << "###################### Parameters ######################" << std::endl;
  return output.str();
}

}

}  // namespace toppic
