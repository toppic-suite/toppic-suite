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
#include <QScrollBar>

#include "common/util/file_util.hpp"
#include "common/util/version.hpp"
#include "common/util/mem_check.hpp"

#include "console/topindex_argument.hpp"

#include "gui/util/command.hpp"
#include "gui/util/gui_message.hpp"

#include "gui/topindex/ui_topindexdialog.h"
#include "gui/topindex/topindexdialog.hpp"


TopIndexDialog::TopIndexDialog(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::TopIndexDialog) {
      arguments_ = toppic::TopIndexArgument::initArguments();
      ui->setupUi(this);
      std::string title = "TopIndex v." + toppic::Version::getVersion();
      QString qstr = QString::fromStdString(title);
      this->setWindowTitle(qstr);
      lastDir_ = ".";
      QFont font;
      QFont outputFont;
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
      font.setFamily(QStringLiteral("Calibri"));
      outputFont.setFamily(QStringLiteral("Consolas"));
#else
      font.setFamily(QStringLiteral("Monospace"));
      outputFont.setFamily(QStringLiteral("Monospace"));
#endif
      font.setPixelSize(12);
      outputFont.setPixelSize(12);
      QApplication::setFont(font);
      ui->outputTextBrowser->setFont(outputFont);

      TopIndexDialog::on_defaultButton_clicked();
    }

TopIndexDialog::~TopIndexDialog() {
  if(process_.state()!=QProcess::NotRunning) {
    process_.kill();
  }
  delete ui;
}

void TopIndexDialog::on_databaseFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a protein database file",
      lastDir_,
      "Database files (*.fasta *.fa)");
  updatedir(s);
  ui->databaseFileEdit->setText(s);
}
void TopIndexDialog::on_fixedModFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a fixed modification file",
      lastDir_,
      "Modification files (*.txt);;All files (*.*)");
  updatedir(s);
  ui->fixedModFileEdit->setText(s);
}

void TopIndexDialog::closeEvent(QCloseEvent *event) {
  if(process_.state()!=QProcess::NotRunning) {
    if (!continueToClose()) {
      event->ignore();
      return;
    }
    else {
      process_.kill();
    }
  }
  event->accept();
  return;
}

void TopIndexDialog::on_clearButton_clicked() {
  ui->databaseFileEdit->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the data.");
  ui->outputButton->setEnabled(false);
}

void TopIndexDialog::on_defaultButton_clicked() {
  arguments_ = toppic::TopIndexArgument::initArguments();
  ui->fixedModComboBox->setCurrentIndex(0);
  ui->errorToleranceEdit_2->setText(QString::fromStdString(arguments_["massErrorTolerance"])); 
  ui->threadNumberEdit->setText(QString::fromStdString(arguments_["threadNumber"]));
  ui->fixedModComboBox->setCurrentIndex(0);
  on_fixedModComboBox_currentIndexChanged(0);
  ui->NONECheckBox->setChecked(true);
  ui->NMECheckBox->setChecked(true);
  ui->NMEACCheckBox->setChecked(true);
  ui->MACCheckBox->setChecked(true);
  ui->decoyCheckBox->setChecked(false);
  ui->outputTextBrowser->setText("Click the Start button to process the data.");
}

void TopIndexDialog::updatedir(QString s) {
  if (!s.isEmpty()) {
    //lastDir_ = s;
    lastDir_ = "";
  }
}
void TopIndexDialog::on_startButton_clicked() {
  lockDialog();
  
  std::map<std::string, std::string> argument = this->getArguments();

  std::string cmd = toppic::command::geneTopIndexCommand(argument); 
  QString q_cmd = QString::fromStdString(cmd);
  q_cmd = q_cmd.trimmed();
  QStringList cmd_list = q_cmd.split(" ");
  QString prog = cmd_list[0];
  cmd_list.removeFirst();

  //qDebug() << "start process ";
  process_.start(prog, cmd_list);
  process_.waitForStarted();
  //qDebug() << "start process finished";

  toppic::GuiMessage guiMsg;
  bool finish = false;
  while (!finish) {
    if(process_.state()==QProcess::NotRunning) {
      finish = true;
    }
    bool ready = process_.waitForReadyRead(100);
    if (ready || finish) {
      //qDebug() << "read finished";
      QByteArray byteArray = process_.readAllStandardOutput();
      QString str = QString(byteArray);
      std::string msg = guiMsg.getMsg(str.toStdString());
      if (msg != "") {
        updateMsg(msg); 
      }
    }
    if (finish) {
      QByteArray byteArray = process_.readAllStandardError();
      QString str = QString(byteArray);
      if (process_.exitStatus() != QProcess::NormalExit) {
        str = str + "\nERROR Quit status: Crashed. \n";
        str = str + "ERROR Quit code: " + QString::number(process_.exitCode()) + ".\n";
      }
      std::string msg = guiMsg.getMsg(str.toStdString());
      if (msg != "") {
        updateMsg(msg); 
      }
      //qDebug() << "Status: " << process_.exitStatus();
      //qDebug() << "Code: " << process_.exitCode();
    }
    sleep(100);
  }
  unlockDialog();
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
  std::string db_file_name = ui->databaseFileEdit->text().toStdString();

  std::string dir = toppic::file_util::directory(db_file_name);
  QString outPath = dir.c_str();

  QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
}

std::map<std::string, std::string> TopIndexDialog::getArguments() {
  QString path = QCoreApplication::applicationFilePath();
  arguments_["executiveDir"] = toppic::file_util::getExecutiveDir(path.toStdString());
  if (toppic::file_util::checkSpace(arguments_["executiveDir"])) {
    ui->outputTextBrowser->setText("Current directory " + QString::fromStdString(arguments_["executiveDir"]) + " contains space and will cause errors in the program!");
  }
  arguments_["resourceDir"] = toppic::file_util::getResourceDir(arguments_["executiveDir"]);
  arguments_["oriDatabaseFileName"] = ui->databaseFileEdit->text().toStdString();

  if (ui->decoyCheckBox->isChecked()) {
    arguments_["searchType"] = "TARGET+DECOY";
    arguments_["databaseFileName"] = arguments_["oriDatabaseFileName"] + "_target_decoy";
  } 
  else {
    arguments_["searchType"] = "TARGET";
    arguments_["databaseFileName"] = arguments_["oriDatabaseFileName"] + "_target";
  }
  arguments_["fixedMod"] = ui->fixedModComboBox->currentText().toStdString();
  if (arguments_["fixedMod"] == "NONE") {
    arguments_["fixedMod"] = "";
  }
  else if (arguments_["fixedMod"] == "Carbamidomethylation on cysteine") {
    arguments_["fixedMod"] = "C57";
  }
  else if (arguments_["fixedMod"] == "Carboxymethylation on cysteine") {
    arguments_["fixedMod"] = "C58";
  }
  if (ui->fixedModComboBox->currentIndex() == 3) {
    arguments_["fixedMod"] = ui->fixedModFileEdit->text().toStdString();
  }
  arguments_["massErrorTolerance"] = ui->errorToleranceEdit_2->text().toStdString();

  arguments_["allowProtMod"] = "";
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
  if (ui->threadNumberEdit->text().toInt() > toppic::mem_check::getMaxThreads("topindex")) {
    int max_thread = toppic::mem_check::getMaxThreads("topindex");
    QMessageBox::StandardButton reply = QMessageBox::warning(this, tr("Warning"),
                         QString("Thread number is too large! Based on the memory size, up to %1 threads can run on this computer. Are you sure you want to proceed?").arg(max_thread).arg(max_thread),
                         QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::No) {
      return true;
    }
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
  QString showInfo = msg.c_str();
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  int vertical_bar_pos = ui->outputTextBrowser->verticalScrollBar()->value();
  int max_bar_pos = ui->outputTextBrowser->verticalScrollBar()->maximum();
  ui->outputTextBrowser->setText(showInfo);
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
  if (max_bar_pos - vertical_bar_pos < 10) {
    vertical_bar_pos = ui->outputTextBrowser->verticalScrollBar()->maximum();
  }
  ui->outputTextBrowser->verticalScrollBar()->setValue(vertical_bar_pos);
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
