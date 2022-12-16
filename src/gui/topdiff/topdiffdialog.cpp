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

#include <sstream>

#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QCloseEvent>
#include <QDesktopServices>
#include <QScrollBar>

#include "common/util/file_util.hpp"
#include "common/base/base_data.hpp"
#include "common/util/version.hpp"

#include "console/topdiff_argument.hpp"

#include "gui/util/command.hpp"
#include "gui/util/gui_message.hpp"

#include "gui/topdiff/ui_topdiffdialog.h"
#include "gui/topdiff/topdiffdialog.hpp"

TopDiffDialog::TopDiffDialog(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::TopDiffDialog) {
      ui->setupUi(this);
      std::string title = "TopDiff v." + toppic::Version::getVersion();
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

      TopDiffDialog::on_defaultButton_clicked();
    }

TopDiffDialog::~TopDiffDialog() {
  if(process_.state()!=QProcess::NotRunning) {
    process_.kill();
  }
  delete ui;
}

void TopDiffDialog::closeEvent(QCloseEvent *event) {
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

void TopDiffDialog::on_clearButton_clicked() {
  ui->listWidget->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the data.");
  ui->outputButton->setEnabled(false);
}

void TopDiffDialog::on_defaultButton_clicked() {
  arguments_ = toppic::TopDiffArgument::initArguments();
  ui->toolComboBox->setCurrentIndex(0);
  ui->precErrorEdit->setText(QString::fromStdString(arguments_["errorTolerance"]));
  ui->outputEdit->setText(QString::fromStdString(arguments_["mergedOutputFileName"])); 
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
      "Spectrum files (*ms2.msalign)");
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
    //lastDir_ = s;
    lastDir_ = "";
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
  lockDialog();
  std::map<std::string, std::string> argument = this->getArguments();
  std::vector<std::string> spec_file_lst = this->getSpecFileList();

  std::string cmd = toppic::command::geneTopDiffCommand(argument, spec_file_lst_);
  QString q_cmd = QString::fromStdString(cmd);
  q_cmd = q_cmd.trimmed();
  QStringList cmd_list = q_cmd.split(" ");
  QString prog = cmd_list[0];
  cmd_list.removeFirst();

  process_.start(prog, cmd_list);
  process_.waitForStarted();

  toppic::GuiMessage guiMsg;
  bool finish = false;
  while (!finish) {
    if(process_.state()==QProcess::NotRunning) {
      finish = true;
    }
    bool ready = process_.waitForReadyRead(100);
    if (ready || finish) {
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
  if (toppic::file_util::checkSpace(exe_dir)) {
    ui->outputTextBrowser->setText("Current directory " + QString::fromStdString(exe_dir) + " contains space and will cause errors in the program!");
  }
  arguments_["resourceDir"] = toppic::file_util::getResourceDir(exe_dir);
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
  
  ui->outputEdit->setEnabled(false);
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

  ui->precErrorEdit->setEnabled(true);
  ui->outputEdit->setEnabled(true);
  ui->toolComboBox->setEnabled(true);
}

bool TopDiffDialog::checkError() {
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

void TopDiffDialog::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait) {
    QCoreApplication::processEvents();
  }
}

