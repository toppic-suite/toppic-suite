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

#ifndef TOPPIC_TOPFD_FEATURE_DETECT_FEATURE_DETECT_HPP_
#define TOPPIC_TOPFD_FEATURE_DETECT_FEATURE_DETECT_HPP_

#include <string>

namespace toppic {

namespace feature_detect {

void process(int frac_id, const std::string &sp_file_name, 
             bool miss_level_one, const std::string &resource_dir);
 
};

}

#endif
