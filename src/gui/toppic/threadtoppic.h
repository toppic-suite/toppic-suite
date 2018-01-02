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

#ifndef PROT_GUI_THREADTOPPIC_H
#define PROT_GUI_THREADTOPPIC_H

#include <map>
#include <string>

#include <QThread>

namespace Ui {
class threadtoppic;
}

class threadtoppic : public QThread {
 Q_OBJECT
 public:
  explicit threadtoppic(QObject* par);
  ~threadtoppic() {}
  void run();
  void setPar(std::map<std::string, std::string> arguments) {
    arguments_ = arguments;
  }
 private:
  std::map<std::string, std::string> arguments_;
};

#endif  // PROT_GUI_THREADTOPPIC_H
