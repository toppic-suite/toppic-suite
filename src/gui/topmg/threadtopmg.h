//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#ifndef PROT_GUI_THREADTOPMG_H
#define PROT_GUI_THREADTOPMG_H

#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <ctime>

#include "common/util/version.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"

#include "spec/msalign_util.hpp"
#include "spec/feature_util.hpp"

#include "prsm/prsm_util.hpp"

#include "console/topmg_argument.hpp"
#include "console/topmg_process.hpp"

#include <QThread>

namespace Ui {
class threadtopmg;
}

class threadtopmg : public QThread {
 Q_OBJECT
 public:
  explicit threadtopmg(QObject* par) : QThread(par) {}

  ~threadtopmg() {}

  void run();

  void setPar(std::map<std::string, std::string> arguments,
             const std::vector<std::string> & spec_file_lst) {
    arguments_ = arguments;
    spec_file_lst_ = spec_file_lst;
  }

 private:
  std::map<std::string, std::string> arguments_;

  std::vector<std::string> spec_file_lst_;
};

#endif  // PROT_GUI_THREADTOPMG_H
