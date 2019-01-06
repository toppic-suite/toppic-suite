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


#ifndef TOPPIC_GUI_THREADTOPMERGE_H
#define TOPPIC_GUI_THREADTOPMERGE_H

#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

#include <QThread>

namespace Ui {
class ThreadTopMerge;
}

class ThreadTopMerge : public QThread {
  Q_OBJECT

 public:
  explicit ThreadTopMerge(QObject* par) : QThread(par) {}

  ~ThreadTopMerge() {}

  void run();

  void setPar(std::map<std::string, std::string> arguments,
              const std::vector<std::string> & proteoform_file_lst) {
    arguments_ = arguments;
    proteoform_file_lst_ = proteoform_file_lst;
  }

 private:
  std::map<std::string, std::string> arguments_;

  std::vector<std::string> proteoform_file_lst_;
};

#endif

