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


#ifndef PROT_GUI_THREADTOPFD_H
#define PROT_GUI_THREADTOPFD_H

#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

#include <QThread>

#include "util/file_util.hpp"
#include "util/string_util.hpp"
#include "console/topfd_argument.hpp"

#include "feature/topfd_process.hpp"

namespace Ui {
class ThreadTopFD;
}

class ThreadTopFD : public QThread {
  Q_OBJECT

 public:
  explicit ThreadTopFD(QObject* par) : QThread(par) {}

  ~ThreadTopFD() {}

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

#endif

