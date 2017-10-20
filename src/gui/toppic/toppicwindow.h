// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef PROT_GUI_TOPPICWINDOW_H
#define PROT_GUI_TOPPICWINDOW_H

#include <QMainWindow>
#include <QMouseEvent>
#include "threadtoppic.h"
#include <map>
#include <string>

namespace Ui {
class toppicWindow;
}

class toppicWindow : public QMainWindow {
 Q_OBJECT

 public:
  explicit toppicWindow(QWidget *parent = 0);
  ~toppicWindow();

  private slots:
      void on_databaseFileButton_clicked();
  void on_spectrumFileButton_clicked();
  void on_fixedModFileButton_clicked();
  void on_modFileButton_clicked();
  void on_topfdFeatureFileButton_clicked();
  void on_clearButton_clicked();
  void on_defaultButton_clicked();
  void on_startButton_clicked();
  void on_exitButton_clicked();
  void on_outputButton_clicked();
  void on_fixedModComboBox_currentIndexChanged(int index);
  void on_topfdFeatureCheckBox_clicked(bool checked);
  void on_errorToleranceEdit_textChanged(QString string);
  void on_generatingFunctionCheckBox_clicked(bool checked);
  void on_NONECheckBox_clicked(bool checked);
  void on_NMECheckBox_clicked(bool checked);
  void on_NMEACCheckBox_clicked(bool checked);
  void on_MACCheckBox_clicked(bool checked);
  void on_numModComboBox_currentIndexChanged(int index);
  void on_cutoffSpectralTypeComboBox_currentIndexChanged(int index);
  void on_cutoffProteoformTypeComboBox_currentIndexChanged(int index);
  void on_decoyCheckBox_clicked(bool checked);

 private:
  Ui::toppicWindow *ui;
  QString lastDir_;
  std::map<std::string, std::string> arguments_;
  void initArguments();
  std::map<std::string, std::string> getArguments();
  void lockDialog();
  void unlockDialog();
  bool checkError();
  void updateMsg(std::string msg);
  void updatedir(QString s);
  void showArguments();
  void sleep(int wait);
  threadtoppic* thread_;
  QString showInfo;
  void closeEvent(QCloseEvent *event);
  bool continueToClose();
  bool nterminalerror();
  bool event(QEvent *event);
};

#endif  // PROT_GUI_TOPPICWINDOW_H
