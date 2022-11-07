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

#ifndef TOPPIC_GUI_TOPINDEX_DIALOG_HPP
#define TOPPIC_GUI_TOPINDEX_DIALOG_HPP

#include <map>
#include <string>

#include <QMainWindow>
#include <QProcess>

namespace Ui {
class TopIndexDialog;
}

class TopIndexDialog : public QMainWindow {
  Q_OBJECT

public:
  explicit TopIndexDialog(QWidget *parent = 0);
  ~TopIndexDialog();

private slots:
  void on_clearButton_clicked();

  void on_defaultButton_clicked();

  void on_startButton_clicked();

  void on_exitButton_clicked();

  void on_outputButton_clicked();

  void on_databaseFileButton_clicked();

  void on_fixedModFileButton_clicked();

  void on_fixedModComboBox_currentIndexChanged(int index);

   void on_NONECheckBox_clicked(bool checked);

  void on_NMECheckBox_clicked(bool checked);

  void on_NMEACCheckBox_clicked(bool checked);

  void on_MACCheckBox_clicked(bool checked);
  
  bool nterminalerror();

private:
  QString lastDir_;

  std::map<std::string, std::string> arguments_;

  std::vector<std::string> spec_file_lst_;

  QProcess process_;

  Ui::TopIndexDialog *ui;

  std::map<std::string, std::string> getArguments();

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
