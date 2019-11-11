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

#ifndef TOPPIC_GUI_TOPFD_THREADTOPFD_HPP
#define TOPPIC_GUI_TOPFD_THREADTOPFD_HPP

#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

#include <QThread>

#include "topfd/common/topfd_para.hpp"

namespace Ui {
class ThreadTopFD;
}

class ThreadTopFD : public QThread {
  Q_OBJECT

 public:
  explicit ThreadTopFD(QObject* par) : QThread(par) {}

  ~ThreadTopFD() {}

  void run();

  void setPar(toppic::TopfdParaPtr para_ptr, 
              const std::vector<std::string> & spec_file_lst); 

 private:
  toppic::TopfdParaPtr para_ptr_;

  std::vector<std::string> spec_file_lst_;
};

#endif

