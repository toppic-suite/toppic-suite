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

#ifndef GUI_UTIL_GUI_MESSAGE_HPP
#define GUI_UTIL_GUI_MESSAGE_HPP

#include <string>
#include <sstream>

#include <QProcess>
#include <QTextBrowser>

namespace toppic {

class GuiMessage {
 public:
  GuiMessage();
  std::string getMsg(std::string new_msg);

 private:
  std::stringstream buffer_;
  std::string info_;
  int processed_len_ = 0;
  std::string processed_lines_ = ""; 
  std::string current_line_ = "";
  unsigned cursor_pos_ = 0;
  bool finish_ = false;
};

}
#endif
