//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "gui/topdiff/topdiffdialog.h"
#include "gui/topdiff/ui_topdiffdialog.h"
#include "gui/topdiff/threadtopdiff.h"


TopDiffDialog::TopDiffDialog(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::TopDiffDialog) {
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
      thread_ = new ThreadTopDiff(this);
      showInfo = "";
      TopDiffDialog::on_defaultButton_clicked();
    }

TopDiffDialog::~TopDiffDialog() {
  thread_->terminate();
  thread_->wait();
  delete ui;
}

void TopDiffDialog::on_databaseFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a protein database file",
      lastDir_,
      "Database files(*.fasta *.fa)");
  updatedir(s);
  ui->databaseFileEdit->setText(s);
}

void TopDiffDialog::on_fixedModFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a fixed modification file",
      lastDir_,
      "Modification files(*.txt);;All files(*.*)");
  updatedir(s);
  ui->fixedModFileEdit->setText(s);
}

void TopDiffDialog::on_fixedModComboBox_currentIndexChanged(int index) {
  if (index == 3) {
    ui->fixedModFileEdit->setEnabled(true);
    ui->fixedModFileButton->setEnabled(true);
  } else {
    ui->fixedModFileEdit->setEnabled(false);
    ui->fixedModFileButton->setEnabled(false);
  }
}


void TopDiffDialog::closeEvent(QCloseEvent *event) {
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

void TopDiffDialog::initArguments() {
  arguments_["executiveDir"] = "";
  arguments_["resourceDir"] = "";
  arguments_["databaseFileName"] = "";
  arguments_["fixedMod"] = "";
  arguments_["errorTolerance"] = "1.2";
  arguments_["toolName"] = "toppic";
  arguments_["mergedOutputFileName"] = "sample_diff.csv";
}

void TopDiffDialog::on_clearButton_clicked() {
  ui->databaseFileEdit->clear();
  ui->listWidget->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the data.");
  ui->outputButton->setEnabled(false);
}

void TopDiffDialog::on_defaultButton_clicked() {
  ui->fixedModFileEdit->clear();
  ui->fixedModComboBox->setCurrentIndex(0);
  on_fixedModComboBox_currentIndexChanged(0);
  ui->toolComboBox->setCurrentIndex(0);
  ui->precErrorEdit->setText("1.2");
  ui->outputEdit->setText("sample_diff.csv");
  ui->outputTextBrowser->setText("Click the Start button to process the data.");
}

std::vector<std::string> TopDiffDialog::getSpecFileList() {
  spec_file_lst_.clear();
  for (int i = 0; i < ui->listWidget->count(); i++) {
    spec_file_lst_.push_back(ui->listWidget->item(i)->text().toStdString());
  }
  return spec_file_lst_;
}

void TopDiffDialog::on_addButton_clicked() {
  QStringList idfiles = QFileDialog::getOpenFileNames(
      this,
      "Select spectrum files",
      lastDir_,
      "Spectrum files(*ms2.msalign)");
  for (int i = 0; i < idfiles.size(); i++) {
    QString idfile = idfiles.at(i);
    updatedir(idfile);
    if (ableToAdd(idfile)) {
      ui->listWidget->addItem(new QListWidgetItem(idfile));
    }
  }
}

void TopDiffDialog::updatedir(QString s) {
  if (!s.isEmpty()) {
    lastDir_ = s;
  }
}
bool TopDiffDialog::ableToAdd(QString idfile) {
  bool able = true;
  if (idfile != "") {
    if (idfile.toStdString().length() > 200) {
      QMessageBox::warning(this, tr("Warning"),
                           tr("The file path is too long!"),
                           QMessageBox::Yes);
      able = false;
    } else {
      for (int i = 0; i < ui->listWidget->count(); i++) {
        if (idfile == ui->listWidget->item(i)->text()) {
          able = false;
        }
      }
    }
  } else {
    able = false;
  }
  return able;
}

void TopDiffDialog::on_delButton_clicked() {
  QListWidgetItem *delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
}

void TopDiffDialog::on_startButton_clicked() {
  std::stringstream buffer;
  std::streambuf *oldbuf = std::cout.rdbuf(buffer.rdbuf());
  if (checkError()) {
    return;
  }
  lockDialog();
  ui->outputTextBrowser->setText(showInfo);
  std::map<std::string, std::string> argument = this->getArguments();
  std::vector<std::string> spec_file_lst = this->getSpecFileList();
  thread_->setPar(argument, spec_file_lst);
  thread_->start();

  std::string info;
  int processed_len = 0;
  std::string processed_lines = ""; 
  std::string current_line = "";
  unsigned cursor_pos = 0;

  while (true) {
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
    if (thread_->isFinished()) {
      break;
    }
    sleep(100);
  }
  unlockDialog();
  showInfo = "";
  thread_->exit();
  std::cout.rdbuf(oldbuf);
}

void TopDiffDialog::on_exitButton_clicked() {
  close();
}

bool TopDiffDialog::continueToClose() {
  if (QMessageBox::question(this,
                            tr("Quit"),
                            tr("TopDiff is still running. Are you sure you want to quit?"),
                            QMessageBox::Yes | QMessageBox::No,
                            QMessageBox::No)
      == QMessageBox::Yes) {
    return true;
  } else {
    return false;
  }
}

void TopDiffDialog::on_outputButton_clicked() {
  std::string spec_file_name = "";
  if (spec_file_lst_.size() > 0) {
    spec_file_name = spec_file_lst_[0];
  }
  std::string dir = toppic::file_util::directory(spec_file_name);
  QString outPath = dir.c_str();
  QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
}

std::map<std::string, std::string> TopDiffDialog::getArguments() {
  QString path = QCoreApplication::applicationFilePath();
  std::string exe_dir = toppic::file_util::getExecutiveDir(path.toStdString());
  arguments_["executiveDir"] = exe_dir;
  arguments_["resourceDir"] = toppic::file_util::getResourceDir(exe_dir);
  arguments_["databaseFileName"] = ui->databaseFileEdit->text().toStdString();

  arguments_["fixedMod"] = ui->fixedModComboBox->currentText().toStdString();
  if (arguments_["fixedMod"] == "NONE") {
    arguments_["fixedMod"] = "";
  }
  if (ui->fixedModComboBox->currentIndex() == 3) {
    arguments_["fixedMod"] = ui->fixedModFileEdit->text().toStdString();
  }
  arguments_["toolName"] = ui->toolComboBox->currentText().toStdString();
  arguments_["errorTolerance"] = ui->precErrorEdit->text().toStdString();
  arguments_["mergedOutputFileName"] = ui->outputEdit->text().trimmed().toStdString();
  return arguments_;
}

void TopDiffDialog::lockDialog() {
  ui->addButton->setEnabled(false);
  ui->delButton->setEnabled(false);
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->outputButton->setEnabled(false);
  
  ui->databaseFileButton->setEnabled(false);
  ui->databaseFileEdit->setEnabled(false);
  ui->outputEdit->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  ui->fixedModFileEdit->setEnabled(false);
  ui->fixedModFileButton->setEnabled(false);
  ui->toolComboBox->setEnabled(false);
  ui->precErrorEdit->setEnabled(false);
}

void TopDiffDialog::unlockDialog() {
  ui->addButton->setEnabled(true);
  ui->delButton->setEnabled(true);
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);

  ui->databaseFileButton->setEnabled(true);
  ui->databaseFileEdit->setEnabled(true);
  ui->fixedModFileEdit->setEnabled(true);
  ui->precErrorEdit->setEnabled(true);
  ui->outputEdit->setEnabled(true);
  ui->fixedModComboBox->setEnabled(true);
  on_fixedModComboBox_currentIndexChanged(ui->fixedModComboBox->currentIndex());
  ui->toolComboBox->setEnabled(true);
}

bool TopDiffDialog::checkError() {
  if (ui->databaseFileEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Database file is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->listWidget->count() == 0) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Spectrum files are not selected!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->precErrorEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Error tolerance is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->outputEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Output filename is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  return false;
}

void TopDiffDialog::updateMsg(std::string msg) {
  showInfo = msg.c_str();
  ui->outputTextBrowser->setText(showInfo);
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
}

void TopDiffDialog::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait) {
    QCoreApplication::processEvents();
  }
}

