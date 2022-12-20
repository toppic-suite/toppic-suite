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
#include <iostream>

#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QCloseEvent>
#include <QDesktopServices>
#include <QScrollBar>
#include <QProcess>
#include <QDebug>

#include "common/util/file_util.hpp"
#include "common/util/version.hpp"
#include "common/util/mem_check.hpp"

#include "gui/util/command.hpp"
#include "gui/util/gui_message.hpp"
#include "gui/topfd/ui_topfddialog.h"
#include "gui/topfd/topfddialog.hpp"

TopFDDialog::TopFDDialog(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::TopFDDialog) {
      para_ptr_ = std::make_shared<toppic::TopfdPara>();
      ui->setupUi(this);
      std::string title = "TopFD v." + toppic::Version::getVersion();
      QString qstr = QString::fromStdString(title);
      this->setWindowTitle(qstr);
      lastDir_ = ".";
      ui->maxChargeEdit->setValidator(new QIntValidator(1, 100, this));
      ui->maxMassEdit->setValidator(new QIntValidator(1, 1000000, this));
      QRegExp rx1("^0\\.[0]\\d{0,2}[1-9]|0.1$");
      QRegExpValidator *validator1 = new QRegExpValidator(rx1, this);
      ui->mzErrorEdit->setValidator(validator1);
      QRegExp rx2("^\\d{1,6}\\.\\d{0,2}$");
      QRegExpValidator *validator2 = new QRegExpValidator(rx2, this);
      ui->ms1snRatioEdit->setValidator(validator2);
      ui->ms2snRatioEdit->setValidator(validator2);
      QRegExp rx3("^\\d{1,4}\\.\\d{0,2}|10000$");
      ui->threadNumberEdit->setValidator(new QIntValidator(0, 2147483647, this));
      QRegExpValidator *validator3 = new QRegExpValidator(rx3, this);
      ui->windowSizeEdit->setValidator(validator3);
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
      TopFDDialog::on_defaultButton_clicked();
    }

TopFDDialog::~TopFDDialog() {
  if(process_.state()!=QProcess::NotRunning) {
    process_.kill();
  }
  delete ui;
}

void TopFDDialog::closeEvent(QCloseEvent *event) {
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

void TopFDDialog::on_clearButton_clicked() {
  ui->listWidget->clear();
  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
  lastDir_ = "/";
}

void TopFDDialog::on_defaultButton_clicked() {
  // default para_ptr
  toppic::TopfdParaPtr para_ptr = std::make_shared<toppic::TopfdPara>();
  ui->maxChargeEdit->setText(QString::number(para_ptr->getMaxCharge()));
  ui->maxMassEdit->setText(QString::number(para_ptr->getMaxMass()));
  ui->mzErrorEdit->setText(QString::number(para_ptr->getMzError()));
  ui->ms1snRatioEdit->setText(QString::number(para_ptr->getMsOneSnRatio()));
  ui->ms2snRatioEdit->setText(QString::number(para_ptr->getMsTwoSnRatio()));
  ui->windowSizeEdit->setText(QString::number(para_ptr->getPrecWindow()));
  ui->threadNumberEdit->setText(QString::number(para_ptr->getThreadNum()));
  ui->envCNNCheckBox->setChecked(para_ptr->isUseEnvCnn());
  ui->missLevelOneCheckBox->setChecked(para_ptr->isMissingLevelOne());
  ui->geneHTMLCheckBox->setChecked(para_ptr->isGeneHtmlFolder());
  ui->disableFilteringCheckBox->setChecked(!para_ptr->isDoFinalFiltering());

  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
  ui->activationComboBox->setCurrentIndex(0);
}

std::vector<std::string> TopFDDialog::getSpecFileList() {
  spec_file_lst_.clear();
  for (int i = 0; i < ui->listWidget->count(); i++) {
    spec_file_lst_.push_back(ui->listWidget->item(i)->text().toStdString());
  }
  return spec_file_lst_;
}

void TopFDDialog::on_addButton_clicked() {
  QStringList spfiles = QFileDialog::getOpenFileNames(
      this,
      "Select spectrum files",
      lastDir_,
      "Spectra files (*.mzXML *.mzML *.mzxml *.mzml)");
  for (int i = 0; i < spfiles.size(); i++) {
    QString spfile = spfiles.at(i);
    updatedir(spfile);
    if (ableToAdd(spfile)) {
      ui->listWidget->addItem(new QListWidgetItem(spfile));
    }
  }
}

void TopFDDialog::updatedir(QString s) {
  if (!s.isEmpty()) {
    //lastDir_ = s;
    lastDir_ = "";
  }
}
bool TopFDDialog::ableToAdd(QString spfile) {
  bool able = true;
  if (spfile != "") {
    if (spfile.toStdString().length() > 200) {
      QMessageBox::warning(this, tr("Warning"),
                           tr("The spectrum file path is too long!"),
                           QMessageBox::Yes);
      able = false;
    } else {
      for (int i = 0; i < ui->listWidget->count(); i++) {
        if (spfile == ui->listWidget->item(i)->text()) {
          able = false;
        }
      }
    }
  } else {
    able = false;
  }
  return able;
}

void TopFDDialog::on_delButton_clicked() {
  QListWidgetItem *delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
}

void TopFDDialog::on_startButton_clicked() {
  lockDialog();
  toppic::TopfdParaPtr paraPtr = this->getParaPtr();
  std::vector<std::string> specFileList = this->getSpecFileList();

  std::string cmd = toppic::command::geneTopfdCommand(para_ptr_, spec_file_lst_);
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

void TopFDDialog::on_exitButton_clicked() {
  close();
}

bool TopFDDialog::continueToClose() {
  if (QMessageBox::question(this,
                            tr("Quit"),
                            tr("TopFD is still running. Are you sure you want to quit?"),
                            QMessageBox::Yes | QMessageBox::No,
                            QMessageBox::No)
      == QMessageBox::Yes) {
    return true;
  } else {
    return false;
  }
}

void TopFDDialog::on_outputButton_clicked() {
  std::string sp_file_name = "";
  if (spec_file_lst_.size() > 0) {
    sp_file_name = spec_file_lst_[0];
  }
  std::string dir = toppic::file_util::directory(sp_file_name);
  QString outPath = dir.c_str();
  QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
}

toppic::TopfdParaPtr TopFDDialog::getParaPtr() {
  QString path = QCoreApplication::applicationFilePath();
  std::string exe_dir = toppic::file_util::getExecutiveDir(path.toStdString());
  para_ptr_->setExeDir(exe_dir);
  if (toppic::file_util::checkSpace(exe_dir)) {
    ui->outputTextBrowser->setText("Current directory " + QString::fromStdString(exe_dir) + " contains space and will cause errors in the program!");
  }
  para_ptr_->setResourceDir(toppic::file_util::getResourceDir(exe_dir));
  para_ptr_->setMaxCharge(std::stoi(ui->maxChargeEdit->text().toStdString()));
  para_ptr_->setMaxMass(std::stod(ui->maxMassEdit->text().toStdString()));
  para_ptr_->setMzError(std::stod(ui->mzErrorEdit->text().toStdString()));
  para_ptr_->setMsOneSnRatio(std::stod(ui->ms1snRatioEdit->text().toStdString()));
  para_ptr_->setMsTwoSnRatio(std::stod(ui->ms2snRatioEdit->text().toStdString()));
  para_ptr_->setPrecWindow(std::stod(ui->windowSizeEdit->text().toStdString()));
  para_ptr_->setMissingLevelOne(ui->missLevelOneCheckBox->isChecked()); 
  para_ptr_->setThreadNum(std::stoi(ui->threadNumberEdit->text().toStdString()));
  para_ptr_->setGeneHtmlFolder(ui->geneHTMLCheckBox->isChecked());
  para_ptr_->setUseEnvCnn(ui->envCNNCheckBox->isChecked());
  para_ptr_->setActivation(ui->activationComboBox->currentText().toStdString());
  para_ptr_->setDoFinalFiltering(!(ui->disableFilteringCheckBox->isChecked()));

  return para_ptr_;
}

void TopFDDialog::lockDialog() {
  ui->addButton->setEnabled(false);
  ui->delButton->setEnabled(false);
  ui->maxChargeEdit->setEnabled(false);
  ui->maxMassEdit->setEnabled(false);
  ui->mzErrorEdit->setEnabled(false);
  ui->ms1snRatioEdit->setEnabled(false);
  ui->ms2snRatioEdit->setEnabled(false);
  ui->threadNumberEdit->setEnabled(false);
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->missLevelOneCheckBox->setEnabled(false);
  ui->windowSizeEdit->setEnabled(false);
  ui->outputButton->setEnabled(false);
  ui->geneHTMLCheckBox->setEnabled(false);
  ui->envCNNCheckBox->setEnabled(false);
  ui->activationComboBox->setEnabled(false);
  ui->disableFilteringCheckBox->setEnabled(false);
}

void TopFDDialog::unlockDialog() {
  ui->addButton->setEnabled(true);
  ui->delButton->setEnabled(true);
  ui->maxChargeEdit->setEnabled(true);
  ui->maxMassEdit->setEnabled(true);
  ui->mzErrorEdit->setEnabled(true);
  ui->ms1snRatioEdit->setEnabled(true);
  ui->ms2snRatioEdit->setEnabled(true);
  ui->threadNumberEdit->setEnabled(true);
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->missLevelOneCheckBox->setEnabled(true);
  ui->windowSizeEdit->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);
  ui->geneHTMLCheckBox->setEnabled(true);
  ui->envCNNCheckBox->setEnabled(true);
  ui->activationComboBox->setEnabled(true);
  ui->disableFilteringCheckBox->setEnabled(true);
}

bool TopFDDialog::checkError() {
  if (ui->maxChargeEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Maximum charge is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->maxMassEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Maximum mass is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->mzErrorEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("M/z error tolerance is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->ms1snRatioEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS1 S/N ratio is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->ms2snRatioEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS2 S/N ratio is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->windowSizeEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Precursor window size is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->threadNumberEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Thread number is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->threadNumberEdit->text().toInt() > toppic::mem_check::getMaxThreads("topfd")) {
    int max_thread = toppic::mem_check::getMaxThreads("topfd");
    QMessageBox::StandardButton reply = QMessageBox::warning(this, tr("Warning"),
                         QString("Thread number is too large! Based on the memory size, up to %1 threads can run on this computer. Are you sure you want to proceed?").arg(max_thread).arg(max_thread),
                         QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::No) {
      return true;
    }
  }

  return false;
}

void TopFDDialog::updateMsg(std::string msg) {
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  int vertical_bar_pos = ui->outputTextBrowser->verticalScrollBar()->value();
  int max_bar_pos = ui->outputTextBrowser->verticalScrollBar()->maximum();
  QString showInfo = msg.c_str();
  ui->outputTextBrowser->setText(showInfo);
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
  if (max_bar_pos - vertical_bar_pos < 10) {
    vertical_bar_pos = ui->outputTextBrowser->verticalScrollBar()->maximum();
  }
  ui->outputTextBrowser->verticalScrollBar()->setValue(vertical_bar_pos);
}

void TopFDDialog::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait) {
    QCoreApplication::processEvents();
  }
}

