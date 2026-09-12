// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_GUI_TOP_CONVERTER_TOPCONVERTERDIALOG_HPP_
#define TOPPIC_GUI_TOP_CONVERTER_TOPCONVERTERDIALOG_HPP_

#include <QMainWindow>
#include <QProcess>
#include <string>

namespace Ui {
class TopConverterDialog;
}

// Front-end of top_sql_converter: picks one mzML/mzXML file and the two grid
// parameters, then runs the command-line tool and shows its output.
class TopConverterDialog : public QMainWindow {
  Q_OBJECT

 public:
  explicit TopConverterDialog(QWidget* parent = nullptr);
  ~TopConverterDialog();

 private slots:
  void on_browseButton_clicked();

  void on_clearButton_clicked();

  void on_defaultButton_clicked();

  void on_startButton_clicked();

  void on_exitButton_clicked();

  void on_outputButton_clicked();

 private:
  Ui::TopConverterDialog* ui;

  QProcess process_;

  QString lastDir_;

  std::string getSpectrumFileName();

  void lockDialog();

  void unlockDialog();

  bool checkError();

  void updateMsg(std::string msg);

  void sleep(int wait);

  void closeEvent(QCloseEvent* event);

  bool continueToClose();
};

#endif
