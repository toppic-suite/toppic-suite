//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_GUI_TOPDIA_TOPDIADIALOG_HPP
#define TOPPIC_GUI_TOPDIA_TOPDIADIALOG_HPP

#include <map>
#include <string>

#include <QMainWindow>
#include <QProcess>

#include "topfd/common/topfd_para.hpp"
#include "topdia/common/topdia_para.hpp"

namespace Ui {
class TopDIADialog;
}

class TopDIADialog : public QMainWindow {
  Q_OBJECT

public:
  explicit TopDIADialog(QWidget *parent = 0);
  ~TopDIADialog();

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

  toppic::TopfdParaPtr topfd_para_ptr_;
  toppic::TopdiaParaPtr topdia_para_ptr_;

  std::vector<std::string> spec_file_lst_;

  QProcess process_;

  Ui::TopDIADialog *ui;

  void getParaPtr();

  std::vector<std::string> getSpecFileList();

  void lockDialog();

  void unlockDialog();

  bool checkError();

  void updateMsg(std::string msg);  

  void updatedir(QString s);

  void sleep(int wait);

  void closeEvent(QCloseEvent *event);

  bool continueToClose();

  bool ableToAdd(QString spfile);
};

#endif
