//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_GUI_TOPFD_TOPFDDIALOG_HPP
#define TOPPIC_GUI_TOPFD_TOPFDDIALOG_HPP

#include <map>
#include <string>

#include <QMainWindow>

#include "topfd/common/topfd_para.hpp"
#include "gui/topfd/threadtopfd.hpp"

namespace Ui {
class TopFDDialog;
}

class TopFDDialog : public QMainWindow {
  Q_OBJECT

public:
  explicit TopFDDialog(QWidget *parent = 0);
  ~TopFDDialog();

private slots:
  void on_clearButton_clicked();

  void on_defaultButton_clicked();

  void on_startButton_clicked();

  void on_exitButton_clicked();

  void on_outputButton_clicked();

  void on_addButton_clicked();

  void on_delButton_clicked();

private:
  QString lastDir_;

  int percentage_;

  //std::map<std::string, std::string> arguments_;
  toppic::TopfdParaPtr para_ptr_;

  std::vector<std::string> spec_file_lst_;

  Ui::TopFDDialog *ui;

  //void initArguments();

  //std::map<std::string, std::string> getArguments();
  toppic::TopfdParaPtr getParaPtr();

  std::vector<std::string> getSpecFileList();

  void lockDialog();

  void unlockDialog();

  bool checkError();

  void updateMsg(std::string msg);  

  void updatedir(QString s);

  void sleep(int wait);

  ThreadTopFD* thread_;

  QString showInfo;

  void closeEvent(QCloseEvent *event);

  bool continueToClose();

  bool ableToAdd(QString spfile);
};

#endif
