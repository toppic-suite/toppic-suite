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

#include <map>
#include <string>

#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QCloseEvent>
#include <QDesktopServices>


#include "base/file_util.hpp"
#include "base/base_data.hpp"
#include "feature/deconv_para.hpp"
#include "feature/deconv_process.hpp"
#include "feature/feature_detect.hpp"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "topfddialog.h"
#include "ui_topfddialog.h"
#include "threadtopfd.h"

TopFDDialog::TopFDDialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::TopFDDialog) {
  initArguments();
  ui->setupUi(this);
  lastDir_ = "/";
  percentage_ = 0;
  ui->maxChargeEdit->setValidator(new QIntValidator(1, 100, this));
  ui->maxMassEdit->setValidator(new QIntValidator(1, 1000000, this));
  QRegExp rx1("^0\\.[0]\\d{0,2}[1-9]|0.1$");
  QRegExpValidator *validator1 = new QRegExpValidator(rx1, this);
  ui->mzErrorEdit->setValidator(validator1);
  QRegExp rx2("^\\d{1,6}\\.\\d{0,2}$");
  QRegExpValidator *validator2 = new QRegExpValidator(rx2, this);
  ui->snRatioEdit->setValidator(validator2);
  QRegExp rx3("^\\d{1,4}\\.\\d{0,2}|10000$");
  QRegExpValidator *validator3 = new QRegExpValidator(rx3, this);
  ui->windowSizeEdit->setValidator(validator3);
  QFont font;
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  font.setFamily(QStringLiteral("Courier New"));
#else
  font.setFamily(QStringLiteral("Monospace"));
#endif
  ui->outputTextBrowser->setFont(font);
  thread_ = new ThreadTopFD(this);
  showInfo = "";
  TopFDDialog::on_defaultButton_clicked();
}

TopFDDialog::~TopFDDialog() {
  thread_->terminate();
  thread_->wait();
  delete ui;
}
void TopFDDialog::closeEvent(QCloseEvent *event) {
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
void TopFDDialog::initArguments() {
  arguments_["executiveDir"] = "";
  arguments_["spectrumFileName"] = "";
  arguments_["inputType"] = "mzXML";
  arguments_["refinePrecMass"] = "true";
  arguments_["missingLevelOne"] = "false";
  arguments_["maxCharge"] = "30";
  arguments_["maxMass"] = "100000";
  arguments_["mzError"] = "0.02";
  arguments_["snRatio"] = "1.0";
  arguments_["keepUnusedPeaks"] = "false";
  arguments_["outMultipleMass"] = "false";
  arguments_["precWindow"] = "3.0";
}

void TopFDDialog::on_clearButton_clicked() {
  ui->fileEdit->clear();
  ui->maxChargeEdit->clear();
  ui->maxMassEdit->clear();
  ui->mzErrorEdit->clear();
  ui->snRatioEdit->clear();
  ui->windowSizeEdit->clear();
  ui->missLevelOneCheckBox->setChecked(false);
  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum file.");
  lastDir_ = "/";
}

void TopFDDialog::on_defaultButton_clicked() {
  ui->maxChargeEdit->setText("30");
  ui->maxMassEdit->setText("100000");
  ui->mzErrorEdit->setText("0.02");
  ui->snRatioEdit->setText("1.0");
  ui->windowSizeEdit->setText("3.0");
  ui->missLevelOneCheckBox->setChecked(false);
  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum file.");
}

void TopFDDialog::on_fileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
                this,
                "Select a spectrum file",
                lastDir_,
                "Spectra files(*.mzXML *.mzML *.mzxml *.mzml);;mzXML files(*.mzXML *.mzxml);;mzML files(*.mzML *.mzml)");
  if (!s.isEmpty()) {
    lastDir_ = s;
  }
  ui->fileEdit->setText(s);
}

void TopFDDialog::on_startButton_clicked() {
  std::stringstream buffer;
  std::streambuf *oldbuf = std::cout.rdbuf(buffer.rdbuf());
  if (checkError()) {
    return;
  }
  lockDialog();
  ui->outputTextBrowser->setText(showInfo);
  // prot::log_level = 2;
  std::map<std::string, std::string> argument = this->getArguments();
  thread_->setPar(argument);
  thread_->start();
  // prot::deconvProcess(argument);

  std::string lastinfo = "";
  std::string nowinfo = "";
  std::string info;
  while (true) {
    // info = getInfo(i);      // Here is the infomation been shown in the infoBox.
    nowinfo = buffer.str();
    info = nowinfo.substr(lastinfo.length());
    info = info.erase(info.find_last_not_of(" \n\r\t") + 1);
    lastinfo = nowinfo;
    if (info != "") {
      updateMsg(info);
    }
    if (thread_->isFinished()) {
      ui->progressBar->setValue(100);
      break;
    }
    sleep(10);
  }
  unlockDialog();
  showInfo = "";
  thread_->exit();
  std::cout.rdbuf(oldbuf);
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

void TopFDDialog::on_outputButton_clicked()
{
  fs::path full_path(arguments_["spectrumFileName"].c_str());
  QString outPath = full_path.remove_filename().string().c_str();
  QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
}
std::map<std::string, std::string> TopFDDialog::getArguments() {
  QString path = QCoreApplication::applicationFilePath();
  std::string exe_dir = prot::FileUtil::getExecutiveDir(path.toStdString());
  arguments_["executiveDir"] = exe_dir;
  arguments_["spectrumFileName"] = ui->fileEdit->text().toStdString();
  arguments_["maxCharge"] = ui->maxChargeEdit->text().toStdString();
  arguments_["maxMass"] = ui->maxMassEdit->text().toStdString();
  arguments_["mzError"] = ui->mzErrorEdit->text().toStdString();
  arguments_["snRatio"] = ui->snRatioEdit->text().toStdString();
  arguments_["precWindow"] = ui->windowSizeEdit->text().toStdString();
  return arguments_;
}

void TopFDDialog::lockDialog() {
  ui->maxChargeEdit->setEnabled(false);
  ui->maxMassEdit->setEnabled(false);
  ui->mzErrorEdit->setEnabled(false);
  ui->snRatioEdit->setEnabled(false);
  ui->fileEdit->setEnabled(false);
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->fileButton->setEnabled(false);
  ui->missLevelOneCheckBox->setEnabled(false);
  ui->windowSizeEdit->setEnabled(false);
  ui->outputButton->setEnabled(false);
}

void TopFDDialog::unlockDialog() {
  ui->maxChargeEdit->setEnabled(true);
  ui->maxMassEdit->setEnabled(true);
  ui->mzErrorEdit->setEnabled(true);
  ui->snRatioEdit->setEnabled(true);
  ui->fileEdit->setEnabled(true);
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->fileButton->setEnabled(true);
  ui->missLevelOneCheckBox->setEnabled(true);
  ui->windowSizeEdit->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);
}

bool TopFDDialog::checkError() {
  if (ui->fileEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Please select a spectrum file!"),
                         QMessageBox::Yes);
    return true;
  }
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
  if (ui->snRatioEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("S/N ratio is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->windowSizeEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Precursor window size is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  return false;
}

QString TopFDDialog::updatePercentage(QString s) {
  QString per = "wait.";
  if (s.left(6) == "Running") {
    percentage_ = 100;
  } else if (s.right(9) == "finished.") {
    per = s.mid((s.size() - 13));
    if (per.at(0) == tr("\t") || per.at(0) == tr(" ")) {
      per = per.mid(1, 1);
    } else {
      per = per.mid(0, 2);
    }
    if (percentage_ < per.toInt()) {
      percentage_ = per.toInt();
    }
  }
  if (percentage_ > -1 && percentage_ < 101) {
    ui->progressBar->setValue(percentage_);
  }
  return per;
}

void TopFDDialog::updateMsg(std::string msg) {
  if (msg.substr(0, 6) == "Running") {msg = "\n" + msg;}
  QString info = msg.c_str();
  if (msg.at(0) == '\r') {
    int lastloc = info.lastIndexOf("\r");
    info = info.right(info.size() - lastloc);
  }
  updatePercentage(info);
  ui->outputTextBrowser->setText(showInfo + info);
  if (msg.at(0) != '\r') {
    showInfo = ui->outputTextBrowser->toPlainText();
  }
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
  if (msg.substr(0, 6) == "Running") {
    unlockDialog();
    percentage_ = 0;
    return;
  }
  // if (percentage_ == 0) {
  //   ui->outputTextBrowser->setText("Prepairing...");
  // }
}

std::string TopFDDialog::getInfo(int i) {
  QString s = QString::number(i);
  sleep(20);
  if (i < 0) {
    return "************************** Parameters *********************";
  } else if (i < 100) {
    return QString("Processing spectrum_" + s + "0...\t" + s + "% finished.").toStdString();
  } else {
    return "TopFD finished.";
  }
}

void TopFDDialog::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait) {
    QCoreApplication::processEvents();
  }
}


