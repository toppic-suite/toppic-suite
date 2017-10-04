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

#include "base/file_util.hpp"
#include "base/base_data.hpp"
#include "feature/deconv_para.hpp"
#include "feature/deconv_process.hpp"
#include "feature/feature_detect.hpp"

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
      ui->maxChargeEdit->setValidator(new QIntValidator(0, 2147483647, this));
      ui->maxMassEdit->setValidator(new QIntValidator(0, 100000, this));
      ui->mzErrorEdit->setValidator(new QDoubleValidator(0.00, 2147483647.00, 2, this));
      ui->snRatioEdit->setValidator(new QDoubleValidator(0.0, 2147483647.0, 1, this));
      thread_ = new ThreadTopFD(this);
    }

TopFDDialog::~TopFDDialog() {
  thread_->terminate();
  thread_->wait();
  delete ui;
}

void TopFDDialog::initArguments() {
  arguments_["executiveDir"] = "";
  arguments_["spectrumFileName"] = "";
  arguments_["inputType"] = "mzXML";
  arguments_["refinePrecMass"]="true";
  arguments_["missingLevelOne"] = "false";
  arguments_["maxCharge"] = "30";
  arguments_["maxMass"] = "100000";
  arguments_["mzError"] = "0.02";
  arguments_["snRatio"] = "1.0";
  arguments_["keepUnusedPeaks"] = "false";
  arguments_["outMultipleMass"] = "false";
  arguments_["precWindow"] = "2.0";
}

void TopFDDialog::on_clearButton_clicked() {
  ui->fileEdit->clear();
  ui->maxChargeEdit->clear();
  ui->maxMassEdit->clear();
  ui->mzErrorEdit->clear();
  ui->snRatioEdit->clear();
  ui->missLevelOneCheckBox->setChecked(false);
  ui->infolabel->clear();
  ui->infolabel->setText("Press start to process the speatra.");
  lastDir_ = "/";
}

void TopFDDialog::on_defaultButton_clicked() {
  ui->maxChargeEdit->setText("30");
  ui->maxMassEdit->setText("100000");
  ui->mzErrorEdit->setText("0.02");
  ui->snRatioEdit->setText("1.0");
  ui->missLevelOneCheckBox->setChecked(false);
  ui->infolabel->clear();
  ui->infolabel->setText("Press start to process the speatra.");
}

void TopFDDialog::on_fileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a spectra file",
      lastDir_,
      "Spectra files(*.mzXML *.mzML);;mzXML files(*.mzXML);;mzML files(*.mzML)");
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
    info = info.erase(info.find_last_not_of(" \n\r\t")+1);
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
  thread_->exit();
  std::cout.rdbuf(oldbuf);
}

void TopFDDialog::on_exitButton_clicked() {
  thread_->terminate();
  thread_->wait();
  close();
}

std::map<std::string, std::string> TopFDDialog::getArguments() {
  QString path;
  QDir dir;
  path = dir.currentPath();
  arguments_["executiveDir"] = path.toStdString();
  arguments_["spectrumFileName"] = ui->fileEdit->text().toStdString();
  arguments_["maxCharge"] = ui->maxChargeEdit->text().toStdString();
  arguments_["maxMass"] = ui->maxMassEdit->text().toStdString();
  arguments_["mzError"] = ui->mzErrorEdit->text().toStdString();
  arguments_["snRatio"] = ui->snRatioEdit->text().toStdString();
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
}

bool TopFDDialog::checkError() {
  if (ui->fileEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Waring"),
                         tr("Please select a spectra file!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->maxChargeEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Waring"),
                         tr("Max charge is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->maxMassEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Waring"),
                         tr("Max mass is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->mzErrorEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Waring"),
                         tr("Mz errror is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->snRatioEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Waring"),
                         tr("S/N ratio is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  return false;
}

QString TopFDDialog::updatePercentage(QString s) {
  QString per = "wait.";
  if (s == "TopFD finished.") {
    percentage_ = 100;
  } else if (s.right(9) == "finished.") {
    per = s.mid((s.size()-13));
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
  QString info = msg.c_str();
  updatePercentage(info);
  ui->infolabel->setText(info);
  if (info == "TopFD finished.") {
    unlockDialog();
    percentage_ = 0;
    return;
  }
  if (percentage_ == 0) {
    ui->infolabel->setText("Prepairing...");
  }
}

std::string TopFDDialog::getInfo(int i) {
  QString s = QString::number(i);
  sleep(20);
  if (i < 0) {
    return "************************** Parameters *********************";
  } else if (i < 100) {
    return QString("Processing spectrum_"+s+"0...\t"+s+"% finished.").toStdString();
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


