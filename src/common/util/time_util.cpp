//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <ctime>

#include "common/util/logger.hpp"
#include "common/util/time_util.hpp"

namespace toppic {

namespace time_util {

std::string getTimeStr() {
  std::time_t cur_time = time(nullptr);
  char buf[50];
  std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&cur_time));
  std::string time_str = buf;
  return time_str;
}

}  // namespace time_util

}  // namespace toppic

