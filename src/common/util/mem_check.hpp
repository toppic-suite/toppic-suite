//Copyright (c) 2014 - 2021, The Trustees of Indiana University.
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

#ifndef TOPPIC_COMMON_UTIL_MEMCHECK_HPP_
#define TOPPIC_COMMON_UTIL_MEMCHECK_HPP_

#include <map>
#include <string>

namespace toppic {
namespace mem_check {

int getMaxThreads(std::string app_name);

double getTotalMemInGb (); 

double getAvailMemInGb ();

bool checkThreadNum(int thread_num, std::string prog); 

}
}
#endif
