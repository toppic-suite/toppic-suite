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

#include <map>
#include <string>
#include <vector>
#include <sstream>

#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QCloseEvent>
#include <QDesktopServices>

#include "common/util/file_util.hpp"
#include "common/base/base_data.hpp"

#include "gui/topindex/topindexdialog.h"
#include "gui/topindex/ui_topindexdialog.h"
#include "gui/topindex/threadtopindex.h"


TopIndexDialog::TopIndexDialog(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::TopIndexDialog) {
      initArguments();
      ui->setupUi(this);
      lastDir_ = ".";
      QFont font;
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
      font.setFamily(QStringLiteral("Calibri"));
#else
      font.setFamily(QStringLiteral("Monospace"));
#endif
      font.setPixelSize(12);
      QApplication::setFont(font);
      ui->outputTextBrowser->setFont(font);
      thread_ = new ThreadTopIndex(this);
      showInfo = "";
      TopIndexDialog::on_defaultButton_clicked();
    }

TopIndexDialog::~TopIndexDialog() {
  thread_->terminate();
  thread_->wait();
  delete ui;
}

void TopIndexDialog::on_databaseFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a protein database file",
      lastDir_,
      "Database files(*.fasta *.fa)");
  updatedir(s);
  ui->databaseFileEdit->setText(s);
}
void TopIndexDialog::on_fixedModFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a fixed modification file",
      lastDir_,
      "Modification files(*.txt);;All files(*.*)");
  updatedir(s);
  ui->fixedModFileEdit->setText(s);
}

void TopIndexDialog::closeEvent(QCloseEvent *event) {
  if (thread_->isRunning()) {
    if (!continueToClose()) {
      event->ignore();
      return;
    }
  }
  thread_->terminate();
  thread_->wait();
  event->accept();
  return;
}

void TopIndexDialog::initArguments() {
  arguments_["oriDatabaseFileName"]="";
  arguments_["databaseFileName"] = "";
  arguments_["databaseBlockSize"] = "1000000";
  arguments_["spectrumFileName"] = "";
  arguments_["activation"] = "FILE";
  arguments_["searchType"] = "TARGET";
  arguments_["fixedMod"] = "";
  arguments_["ptmNumber"] = "1";
  arguments_["massErrorTolerance"] = "15";
  arguments_["allowProtMod"] = "";
  arguments_["numOfTopPrsms"] = "1";
  arguments_["executiveDir"] = ".";
  arguments_["resourceDir"] = "";
  arguments_["filteringResultNumber"] = "20";
  arguments_["threadNumber"] = "1";
  arguments_["groupSpectrumNumber"] = "1";
  arguments_["residueModFileName"] = "";
  arguments_["skipList"] = "";
}

void TopIndexDialog::on_clearButton_clicked() {
  ui->databaseFileEdit->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the data.");
  ui->outputButton->setEnabled(false);
}

void TopIndexDialog::on_defaultButton_clicked() {
  ui->fixedModComboBox->setCurrentIndex(0);
  ui->outputTextBrowser->setText("Click the Start button to process the data.");
  ui->errorToleranceEdit_2->setText("15");
  ui->fixedModComboBox->setCurrentIndex(0);
  on_fixedModComboBox_currentIndexChanged(0);
  ui->activationComboBox->setCurrentIndex(0);
  ui->NONECheckBox->setChecked(true);
  ui->NMECheckBox->setChecked(true);
  ui->NMEACCheckBox->setChecked(true);
  ui->MACCheckBox->setChecked(true);
  ui->decoyCheckBox->setChecked(false);
}

void TopIndexDialog::updatedir(QString s) {
  if (!s.isEmpty()) {
    lastDir_ = s;
  }
}

void TopIndexDialog::on_startButton_clicked() {
  std::stringstream buffer;
  std::streambuf *oldbuf = std::cout.rdbuf(buffer.rdbuf());
  if (checkError()) {
    return;
  }
  lockDialog();
  
  ui->outputTextBrowser->setText(showInfo);
  std::map<std::string, std::string> argument = this->getArguments();
  
  thread_->setPar(argument);
  thread_->start();

  std::string info;
  int processed_len = 0;
  std::string processed_lines = ""; 
  std::string current_line = "";
  unsigned cursor_pos = 0;
  bool finish = false;

  while (true) {
    if (thread_->isFinished()) {
      finish = true;
    }
    // Here is the infomation been shown in the infoBox.
    info = buffer.str();
    std::string new_info = info.substr(processed_len);
    processed_len = info.length();

    if (new_info.size() > 0) {
      for (unsigned i = 0; i < new_info.size(); i++) {
        // new line
        if (new_info.at(i) == '\n') {
          processed_lines = processed_lines + current_line + '\n';
          current_line = "";
          cursor_pos = 0;
        }
        // CF
        if (new_info.at(i) == '\r') {
          cursor_pos = 0;
        }
        // add a new charactor
        if (new_info.at(i) != '\n' && new_info.at(i) != '\r') {
          if (cursor_pos < current_line.length()) {
            current_line[cursor_pos] = new_info.at(i);
          }
          else {
            current_line = current_line + new_info.at(i);
          }
          cursor_pos++;
        }
      }
      updateMsg(processed_lines + current_line);
    }
    if (finish) {
      break;
    }
    sleep(100);
  }
  unlockDialog();
  showInfo = "";
  thread_->exit();
  std::cout.rdbuf(oldbuf);
  
}

void TopIndexDialog::on_exitButton_clicked() {
  close();
}

bool TopIndexDialog::continueToClose() {
  if (QMessageBox::question(this,
                            tr("Quit"),
                            tr("TopIndex is still running. Are you sure you want to quit?"),
                            QMessageBox::Yes | QMessageBox::No,
                            QMessageBox::No)
      == QMessageBox::Yes) {
    return true;
  } else {
    return false;
  }
}

void TopIndexDialog::on_outputButton_clicked() {
  std::string spec_file_name = "";
  if (spec_file_lst_.size() > 0) {
    spec_file_name = spec_file_lst_[0];
  }
  std::string dir = toppic::file_util::directory(spec_file_name);
  QString outPath = dir.c_str();
  QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
}

std::map<std::string, std::string> TopIndexDialog::getArguments() {
  QString path = QCoreApplication::applicationFilePath();
  std::string exe_dir = toppic::file_util::getExecutiveDir(path.toStdString());
  
  arguments_["executiveDir"] = exe_dir;
  arguments_["resourceDir"] = toppic::file_util::getResourceDir(exe_dir);
  arguments_["oriDatabaseFileName"] = ui->databaseFileEdit->text().toStdString();
  arguments_["activation"] = ui->activationComboBox->currentText().toStdString();

  if (ui->decoyCheckBox->isChecked()) {
  arguments_["searchType"] = "TARGET+DECOY";
  arguments_["databaseFileName"] = arguments_["oriDatabaseFileName"] + "_target_decoy";
  } else {
    arguments_["searchType"] = "TARGET";
    arguments_["databaseFileName"] = arguments_["oriDatabaseFileName"] + "_target";
  }
  arguments_["fixedMod"] = ui->fixedModComboBox->currentText().toStdString();
  if (arguments_["fixedMod"] == "NONE") {
    arguments_["fixedMod"] = "";
  }

  if (ui->fixedModComboBox->currentIndex() == 3) {
    arguments_["fixedMod"] = ui->fixedModFileEdit->text().toStdString();
  }
  arguments_["massErrorTolerance"] = ui->errorToleranceEdit_2->text().toStdString();

  if (ui->NONECheckBox->isChecked()) {
    arguments_["allowProtMod"] = arguments_["allowProtMod"] + "NONE";
  }
  if (ui->NMECheckBox->isChecked()) {
    arguments_["allowProtMod"] = arguments_["allowProtMod"] + ",NME";
  }

  if (ui->NMEACCheckBox->isChecked()) {
    arguments_["allowProtMod"] = arguments_["allowProtMod"] + ",NME_ACETYLATION";
  }
  if (ui->MACCheckBox->isChecked()) {
    arguments_["allowProtMod"] = arguments_["allowProtMod"] + ",M_ACETYLATION";
  }
  std::cout << "argument: " << arguments_["allowProtMod"] << std::endl;
    arguments_["threadNumber"] = ui->threadNumberEdit->text().toStdString();
  

  return arguments_;
}

void TopIndexDialog::lockDialog() {
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->outputButton->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  ui->fixedModFileEdit->setEnabled(false);
  ui->fixedModFileButton->setEnabled(false);
  ui->errorToleranceEdit_2->setEnabled(false);
  ui->databaseFileButton->setEnabled(false);
  ui->databaseFileEdit->setEnabled(false);
  ui->threadNumberEdit->setEnabled(false);
  ui->activationComboBox->setEnabled(false);
  ui->NONECheckBox->setEnabled(false);
  ui->NMECheckBox->setEnabled(false);
  ui->NMEACCheckBox->setEnabled(false);
  ui->MACCheckBox->setEnabled(false);
  ui->decoyCheckBox->setEnabled(false);
}

void TopIndexDialog::unlockDialog() {
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);

  ui->databaseFileButton->setEnabled(true);
  ui->databaseFileEdit->setEnabled(true);
  ui->fixedModFileEdit->setEnabled(true);
  ui->errorToleranceEdit_2->setEnabled(true);
  ui->threadNumberEdit->setEnabled(true);
  ui->fixedModComboBox->setEnabled(true);
  on_fixedModComboBox_currentIndexChanged(ui->fixedModComboBox->currentIndex());
  ui->activationComboBox->setEnabled(true);
  ui->NONECheckBox->setEnabled(true);
  ui->NMECheckBox->setEnabled(true);
  ui->NMEACCheckBox->setEnabled(true);
  ui->MACCheckBox->setEnabled(true);
  ui->decoyCheckBox->setEnabled(true);
}

bool TopIndexDialog::checkError() {
  if (ui->databaseFileEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Database file is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->databaseFileEdit->text().toStdString().length() > 200) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("The protein database file path is too long!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->fixedModFileEdit->text().isEmpty() && ui->fixedModComboBox->currentIndex() == 3) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Please select a fixed modification file!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->errorToleranceEdit_2->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Mass error tolerance is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->threadNumberEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Thread number is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  return false;
}

void TopIndexDialog::updateMsg(std::string msg) {
  showInfo = msg.c_str();
  ui->outputTextBrowser->setText(showInfo);
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
}

void TopIndexDialog::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait) {
    QCoreApplication::processEvents();
  }
}
void TopIndexDialog::on_fixedModComboBox_currentIndexChanged(int index) {
  if (index == 3) {
    ui->fixedModFileEdit->setEnabled(true);
    ui->fixedModFileButton->setEnabled(true);
  } else {
    ui->fixedModFileEdit->setEnabled(false);
    ui->fixedModFileButton->setEnabled(false);
  }
}
//void TopIndexDialog::on_errorToleranceEdit_textChanged(QString string) {
  //QString currentText = ui->errorToleranceEdit_2->text();
//}
void TopIndexDialog::on_NONECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NONECheckBox->setChecked(true);
  }
}

void TopIndexDialog::on_NMECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMECheckBox->setChecked(true);
  }
}

void TopIndexDialog::on_NMEACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMEACCheckBox->setChecked(true);
  }
}

void TopIndexDialog::on_MACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->MACCheckBox->setChecked(true);
  }
}
bool TopIndexDialog::nterminalerror() {
  if (ui->NONECheckBox->isChecked() || ui->NMECheckBox->isChecked() || ui->NMEACCheckBox->isChecked() || ui->MACCheckBox->isChecked()) {
    return false;
  } else {
    QMessageBox::warning(this, tr("Warning"),
                         tr("At least one N-terminal form should be selected!"),
                         QMessageBox::Yes);
    return true;
  }
}