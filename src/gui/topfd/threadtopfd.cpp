//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include "base/file_util.hpp"
#include "base/base_data.hpp"
#include "base/logger.hpp"
#include "feature/topfd_process.hpp"

#include "threadtopfd.h"

ThreadTopFD::ThreadTopFD(QObject* par):QThread(par) {}

void ThreadTopFD::run() {
  prot::TopFDProcess(arguments_);
}
