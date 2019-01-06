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

#ifndef TOPPIC_GUI_TOPMERGE_DIALOG_H
#define TOPPIC_GUI_TOPMERGE_DIALOG_H

#include <map>
#include <string>

#include <QMainWindow>

#include "threadtopmerge.h"

namespace Ui {
class TopMergeDialog;
}

class TopMergeDialog : public QMainWindow {
  Q_OBJECT

public:
  explicit TopMergeDialog(QWidget *parent = 0);
  ~TopMergeDialog();

private slots:
  void on_clearButton_clicked();

  void on_defaultButton_clicked();

  void on_startButton_clicked();

  void on_exitButton_clicked();

  void on_outputButton_clicked();

  void on_addButton_clicked();

  void on_delButton_clicked();

  void on_databaseFileButton_clicked();

  void on_fixedModFileButton_clicked();

  void on_fixedModComboBox_currentIndexChanged(int index);

private:
  QString lastDir_;

  int percentage_;

  std::map<std::string, std::string> arguments_;

  std::vector<std::string> proteo_file_lst_;

  Ui::TopMergeDialog *ui;

  void initArguments();

  std::map<std::string, std::string> getArguments();

  std::vector<std::string> getProteoFileList();

  void lockDialog();

  void unlockDialog();

  bool checkError();

  void updateMsg(std::string msg);  

  void updatedir(QString s);

  void sleep(int wait);

  ThreadTopMerge* thread_;

  QString showInfo;

  void closeEvent(QCloseEvent *event);

  bool continueToClose();

  bool ableToAdd(QString spfile);
};

#endif
