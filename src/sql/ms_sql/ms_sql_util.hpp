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

#ifndef TOPPIC_SQL_MS_SQL_MS_SQL_UTIL_HPP_
#define TOPPIC_SQL_MS_SQL_MS_SQL_UTIL_HPP_

#include <string>

namespace toppic {

namespace ms_sql_util {

// The database written for a spectrum file: the file name with its extension
// replaced by _3d.db.
std::string getDbFileName(const std::string& spec_file_name);

// Reads the MS1 peaks of an mzML/mzXML file and writes them to
// getDbFileName(spec_file_name) for the 3D visualization. See MsSqlWriter for
// the database layout and the meaning of mz_size and rt_divider.
void convert(const std::string& spec_file_name, double mz_size,
             double rt_divider);

}  // namespace ms_sql_util

}  // namespace toppic

#endif
