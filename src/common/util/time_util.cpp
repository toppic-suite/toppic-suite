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

#include <ctime>
#include <string>

#include "common/util/time_util.hpp"

namespace toppic {

namespace time_util {

// Returns the current local time, e.g. "Sat Jun 07 15:20:00 2026".
std::string getTimeStr() {
  std::time_t cur_time = std::time(nullptr);

  // std::localtime returns a pointer to a shared static std::tm and is therefore
  // not thread-safe; use the reentrant per-platform variant instead.
  std::tm local_tm;
#if defined(_MSC_VER)
  if (localtime_s(&local_tm, &cur_time) != 0) {
    return std::string();
  }
#else
  if (localtime_r(&cur_time, &local_tm) == nullptr) {
    return std::string();
  }
#endif

  char buf[64];
  // strftime returns 0 with unspecified buffer contents if the result does not
  // fit, so guard against it rather than returning garbage.
  if (std::strftime(buf, sizeof(buf), "%a %b %d %H:%M:%S %Y", &local_tm) == 0) {
    return std::string();
  }
  return std::string(buf);
}

}  // namespace time_util

}  // namespace toppic
