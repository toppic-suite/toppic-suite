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

#ifndef TOPPIC_GUI_TOPMGWINDOW_H
#define TOPPIC_GUI_TOPMGWINDOW_H

#include <QMainWindow>
#include <QMouseEvent>
#include "threadtopmg.h"
#include <map>
#include <string>

namespace Ui {
class topmgWindow;
}

class topmgWindow : public QMainWindow {
 Q_OBJECT

 public:
  explicit topmgWindow(QWidget *parent = 0);

  ~topmgWindow();

 private slots:
  void on_databaseFileButton_clicked();

  void on_fixedModFileButton_clicked();
  
  void on_modFileButton_clicked();
  
  void on_skipListButton_clicked();

  void on_clearButton_clicked();

  void on_defaultButton_clicked();

  void on_startButton_clicked();

  void on_exitButton_clicked();

  void on_outputButton_clicked();

  void on_fixedModComboBox_currentIndexChanged(int index);
    
  void on_NONECheckBox_clicked(bool checked);
    
  void on_NMECheckBox_clicked(bool checked);
    
  void on_NMEACCheckBox_clicked(bool checked);
    
  void on_MACCheckBox_clicked(bool checked);
    
  void on_cutoffSpectralTypeComboBox_currentIndexChanged(int index);
    
  void on_cutoffProteoformTypeComboBox_currentIndexChanged(int index);
    
  void on_decoyCheckBox_clicked(bool checked);
    
  void on_addButton_clicked();
    
  void on_delButton_clicked();

 private:
  Ui::topmgWindow *ui;
  
  QString lastDir_;
  
  std::map<std::string, std::string> arguments_;
  
  std::vector<std::string> spec_file_lst_;
  
  void initArguments();
  
  std::map<std::string, std::string> getArguments();
  
  std::vector<std::string> getSpecFileList();
  
  void lockDialog();
  
  void unlockDialog();
  
  bool checkError();
  
  void updateMsg(std::string msg);
  
  void updatedir(QString s);
  
  void showArguments();
  
  void sleep(int wait);
  
  threadtopmg* thread_;
  
  QString showInfo;
  
  void closeEvent(QCloseEvent *event);
  
  bool continueToClose();
  
  bool nterminalerror();
  
  bool ableToAdd(QString spfile);
};

#endif  // TOPPIC_GUI_TOPMGWINDOW_H
