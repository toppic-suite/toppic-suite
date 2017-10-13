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

#ifndef PROT_GUI_TOPFDDIALOG_H
#define PROT_GUI_TOPFDDIALOG_H

#include <map>
#include <string>

#include <QDialog>

#include "threadtopfd.h"

namespace Ui {
class TopFDDialog;
}

class TopFDDialog : public QDialog {
  Q_OBJECT

 public:
  explicit TopFDDialog(QWidget *parent = 0);
  ~TopFDDialog();

  private slots:
  void on_clearButton_clicked();

  void on_defaultButton_clicked();

  void on_fileButton_clicked();

  void on_startButton_clicked();

  void on_exitButton_clicked();

 private:
  QString lastDir_;
  int percentage_;
  std::map<std::string, std::string> arguments_;
  Ui::TopFDDialog *ui;
  void initArguments();
  std::map<std::string, std::string> getArguments();
  std::string getInfo(int i);
  void lockDialog();
  void unlockDialog();
  bool checkError();
  QString updatePercentage(QString s);
  void updateMsg(std::string msg);
  void sleep(int wait);
  ThreadTopFD* thread_;
  QString showInfo;
  void closeEvent(QCloseEvent *event);
  bool continueToClose();
};

#endif
