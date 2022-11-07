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


#include <QString>

#include "common/util/logger.hpp"
#include "gui/util/gui_message.hpp"


namespace toppic {

GuiMessage::GuiMessage() {}

std::string GuiMessage::getMsg(std::string new_msg) {
  buffer_ << new_msg;
  // Here is the infomation been shown in the infoBox.
  info_ = buffer_.str();
  std::string new_info = info_.substr(processed_len_);
  processed_len_ = info_.length();

  if (new_info.size() > 0) {
    for (unsigned i = 0; i < new_info.size(); i++) {
      // new line
      if (new_info.at(i) == '\n') {
        processed_lines_ = processed_lines_ + current_line_ + '\n';
        current_line_ = "";
        cursor_pos_ = 0;
      }
      // CF
      if (new_info.at(i) == '\r') {
        cursor_pos_ = 0;
      }
      // add a new charactor
      if (new_info.at(i) != '\n' && new_info.at(i) != '\r') {
        if (cursor_pos_ < current_line_.length()) {
          current_line_[cursor_pos_] = new_info.at(i);
        }
        else {
          current_line_ = current_line_ + new_info.at(i);
        }
        cursor_pos_++;
      }
    }
    return (processed_lines_ + current_line_);
  }
  return "";
}

}
