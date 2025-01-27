//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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
#include "console/topdia_argument.hpp"

#include "gui/util/command.hpp"
#include "gui/util/gui_message.hpp"
#include "gui/topdia/ui_topdiadialog.h"
#include "gui/topdia/topdiadialog.hpp"

TopDIADialog::TopDIADialog(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::TopDIADialog) {
  topfd_para_ptr_ = toppic::Argument::getTopfdParaPtrForTopdia();
  topdia_para_ptr_ = std::make_shared<toppic::TopdiaPara>();
  ui->setupUi(this);
  std::string title = "TopDIA v." + toppic::Version::getVersion();
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
  ui->threadNumberEdit->setValidator(new QIntValidator(0, 1000, this));
  ui->ms1MinScanNumEdit->setValidator(new QIntValidator(1, 3, this));
  ui->ms2MinScanNumEdit->setValidator(new QIntValidator(1, 3, this));
  ui->pseudoMinPeakNumEdit->setValidator(new QIntValidator(10, 1000, this));
  ui->ms1EcscoreCutoffEdit->setValidator(
      new QDoubleValidator(0.0, 1.0, 4, this));
  ui->ms2EcscoreCutoffEdit->setValidator(
      new QDoubleValidator(0.0, 1.0, 4, this));
  ui->pseudoScoreCutoffEdit->setValidator(
      new QDoubleValidator(0.0, 1.0, 4, this));
  ui->ms1IntePccCutoffEdit->setValidator(
      new QDoubleValidator(0.0, 1.0, 4, this));
  ui->ms2IntePccCutoffEdit->setValidator(
      new QDoubleValidator(0.0, 1.0, 4, this));
  QRegExp rx3("^\\d{1,4}\\.\\d{0,2}|10000$");
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
      TopDIADialog::on_defaultButton_clicked();
}

TopDIADialog::~TopDIADialog() {
  if(process_.state()!=QProcess::NotRunning) {
    process_.kill();
  }
  delete ui;
}

void TopDIADialog::closeEvent(QCloseEvent *event) {
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

void TopDIADialog::on_clearButton_clicked() {
  ui->listWidget->clear();
  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
  lastDir_ = "/";
}

void TopDIADialog::on_defaultButton_clicked() {
  // default para_ptr
  topfd_para_ptr_ = toppic::Argument::getTopfdParaPtrForTopdia();
  topdia_para_ptr_ = std::make_shared<toppic::TopdiaPara>();
  ui->maxChargeEdit->setText(QString::number(topfd_para_ptr_->getMaxCharge()));
  ui->maxMassEdit->setText(QString::number(topfd_para_ptr_->getMaxMass()));
  ui->mzErrorEdit->setText(QString::number(topfd_para_ptr_->getMzError()));
  ui->ms1snRatioEdit->setText(QString::number(topfd_para_ptr_->getMsOneSnRatio()));
  ui->ms2snRatioEdit->setText(QString::number(topfd_para_ptr_->getMsTwoSnRatio()));
  ui->windowSizeEdit->setText(QString::number(topfd_para_ptr_->getPrecWindowWidth()));
  ui->threadNumberEdit->setText(QString::number(topfd_para_ptr_->getThreadNum()));
  ui->msDeconvCheckBox->setChecked(topfd_para_ptr_->isSortUseMsDeconv());
  ui->geneHTMLCheckBox->setChecked(topfd_para_ptr_->isGeneHtmlFolder());
  ui->disableFilteringCheckBox->setChecked(!topfd_para_ptr_->isAANumBasedFilter());
  ui->singleScanNoiseLevelCheckBox->setChecked(topfd_para_ptr_->isUseSingleScanNoiseLevel());

  //////////////////////////////////////
  ui->ms1EcscoreCutoffEdit->setText(QString::number(topfd_para_ptr_->getMs1EcscoreCutoff()));
  ui->ms2EcscoreCutoffEdit->setText(QString::number(topfd_para_ptr_->getMs2EcscoreCutoff()));
  ui->ms1MinScanNumEdit->setText(QString::number(topfd_para_ptr_->getMs1MinScanNum()));
  ui->ms2MinScanNumEdit->setText(QString::number(topfd_para_ptr_->getMs2MinScanNum()));
  ui->pseudoScoreCutoffEdit->setText(QString::number(topdia_para_ptr_->getPseudoScoreCutoff()));
  ui->pseudoMinPeakNumEdit->setText(QString::number(topdia_para_ptr_->getPseudoMinPeaks()));
  ui->ms1IntePccCutoffEdit->setText(QString::number(topdia_para_ptr_->getMs1SeedEnvInteCorrToleCutoff()));
  ui->ms2IntePccCutoffEdit->setText(QString::number(topdia_para_ptr_->getMs2SeedEnvInteCorrToleCutoff()));
  //////////////////////////////////////

  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
  ui->activationComboBox->setCurrentIndex(0);
}

std::vector<std::string> TopDIADialog::getSpecFileList() {
  spec_file_lst_.clear();
  for (int i = 0; i < ui->listWidget->count(); i++) {
    spec_file_lst_.push_back(ui->listWidget->item(i)->text().toStdString());
  }
  return spec_file_lst_;
}

void TopDIADialog::on_addButton_clicked() {
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

void TopDIADialog::updatedir(QString s) {
  if (!s.isEmpty()) {
    //lastDir_ = s;
    lastDir_ = "";
  }
}
bool TopDIADialog::ableToAdd(QString spfile) {
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

void TopDIADialog::on_delButton_clicked() {
  QListWidgetItem *delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
}

void TopDIADialog::on_startButton_clicked() {
  lockDialog();
  this->getParaPtr();
  std::vector<std::string> specFileList = this->getSpecFileList();

  std::string cmd = toppic::command::geneTopdiaCommand(topfd_para_ptr_, topdia_para_ptr_, spec_file_lst_);
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
      if (process_.exitCode() != 0) {
        str = str + "\nERROR Quit status: Crashed. \n";
        str = str + "ERROR Quit code: " + QString::number(process_.exitCode()) + ".\n";
      }
      std::string msg = guiMsg.getMsg(str.toStdString());
      if (msg != "") {
        updateMsg(msg); 
      }
    }
    sleep(100);
  }
  unlockDialog();
}

void TopDIADialog::on_exitButton_clicked() {
  close();
}

bool TopDIADialog::continueToClose() {
  if (QMessageBox::question(this,
                            tr("Quit"),
                            tr("TopDIA is still running. Are you sure you want to quit?"),
                            QMessageBox::Yes | QMessageBox::No,
                            QMessageBox::No)
      == QMessageBox::Yes) {
    return true;
  } else {
    return false;
  }
}

void TopDIADialog::on_outputButton_clicked() {
  std::string sp_file_name = "";
  if (spec_file_lst_.size() > 0) {
    sp_file_name = spec_file_lst_[0];
  }
  std::string dir = toppic::file_util::directory(sp_file_name);
  QString outPath = dir.c_str();
  QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
}

void TopDIADialog::getParaPtr() {
  QString path = QCoreApplication::applicationFilePath();
  std::string exe_dir = toppic::file_util::getExecutiveDir(path.toStdString());
  topfd_para_ptr_->setExeDir(exe_dir);
  if (toppic::file_util::checkSpace(exe_dir)) {
    ui->outputTextBrowser->setText("Current directory " + QString::fromStdString(exe_dir) + " contains space and will cause errors in the program!");
  }
  topfd_para_ptr_->setResourceDir(toppic::file_util::getResourceDir(exe_dir));
  topfd_para_ptr_->setMaxCharge(std::stoi(ui->maxChargeEdit->text().toStdString()));
  topfd_para_ptr_->setMaxMass(std::stod(ui->maxMassEdit->text().toStdString()));
  topfd_para_ptr_->setMzError(std::stod(ui->mzErrorEdit->text().toStdString()));
  topfd_para_ptr_->setMsOneSnRatio(std::stod(ui->ms1snRatioEdit->text().toStdString()));
  topfd_para_ptr_->setMsTwoSnRatio(std::stod(ui->ms2snRatioEdit->text().toStdString()));
  topfd_para_ptr_->setPrecWindowWidth(std::stod(ui->windowSizeEdit->text().toStdString()));
  topfd_para_ptr_->setThreadNum(std::stoi(ui->threadNumberEdit->text().toStdString()));
  topfd_para_ptr_->setGeneHtmlFolder(ui->geneHTMLCheckBox->isChecked());
  topfd_para_ptr_->setSortUseMsDeconv(ui->msDeconvCheckBox->isChecked());
  topfd_para_ptr_->setActivation(ui->activationComboBox->currentText().toStdString());
  topfd_para_ptr_->setAANumBasedFilter(!(ui->disableFilteringCheckBox->isChecked()));

  //////////////////////////////////////
  topfd_para_ptr_->setMs1EcscoreCutoff(std::stod(ui->ms1EcscoreCutoffEdit->text().toStdString()));
  topfd_para_ptr_->setMs2EcscoreCutoff(std::stod(ui->ms2EcscoreCutoffEdit->text().toStdString()));
  topfd_para_ptr_->setMs1MinScanNum(std::stoi(ui->ms1MinScanNumEdit->text().toStdString()));
  topfd_para_ptr_->setMs2MinScanNum(std::stoi(ui->ms2MinScanNumEdit->text().toStdString()));
  topdia_para_ptr_->setPseudoScoreCutoff(std::stod(ui->pseudoScoreCutoffEdit->text().toStdString()));
  topdia_para_ptr_->setPseudoMinPeaks(std::stoi(ui->pseudoMinPeakNumEdit->text().toStdString()));
  topdia_para_ptr_->setMs1SeedEnvInteCorrToleCutoff(std::stod(ui->ms1IntePccCutoffEdit->text().toStdString()));
  topdia_para_ptr_->setMs2SeedEnvInteCorrToleCutoff(std::stod(ui->ms2IntePccCutoffEdit->text().toStdString()));
  //////////////////////////////////////

  topfd_para_ptr_->setUseSingleScanNoiseLevel(ui->singleScanNoiseLevelCheckBox->isChecked());
}
  
void TopDIADialog::lockDialog() {
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
  ui->windowSizeEdit->setEnabled(false);
  ui->outputButton->setEnabled(false);
  ui->geneHTMLCheckBox->setEnabled(false);
  ui->msDeconvCheckBox->setEnabled(false);
  ui->activationComboBox->setEnabled(false);
  ui->disableFilteringCheckBox->setEnabled(false);
  ui->singleScanNoiseLevelCheckBox->setEnabled(false);

  ui->ms1EcscoreCutoffEdit->setEnabled(false);
  ui->ms2EcscoreCutoffEdit->setEnabled(false);
  ui->ms1MinScanNumEdit->setEnabled(false);
  ui->ms2MinScanNumEdit->setEnabled(false);
  ui->pseudoScoreCutoffEdit->setEnabled(false);
  ui->pseudoMinPeakNumEdit->setEnabled(false);
  ui->ms1IntePccCutoffEdit->setEnabled(false);
  ui->ms2IntePccCutoffEdit->setEnabled(false);
}

void TopDIADialog::unlockDialog() {
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
  ui->windowSizeEdit->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);
  ui->geneHTMLCheckBox->setEnabled(true);
  ui->msDeconvCheckBox->setEnabled(true);
  ui->activationComboBox->setEnabled(true);
  ui->disableFilteringCheckBox->setEnabled(true);
  ui->singleScanNoiseLevelCheckBox->setEnabled(true);

  ui->ms1EcscoreCutoffEdit->setEnabled(true);
  ui->ms2EcscoreCutoffEdit->setEnabled(true);
  ui->ms1MinScanNumEdit->setEnabled(true);
  ui->ms2MinScanNumEdit->setEnabled(true);
  ui->pseudoScoreCutoffEdit->setEnabled(true);
  ui->pseudoMinPeakNumEdit->setEnabled(true);
  ui->ms1IntePccCutoffEdit->setEnabled(true);
  ui->ms2IntePccCutoffEdit->setEnabled(true);

}

bool TopDIADialog::checkError() {
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

  if (ui->ms1MinScanNumEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS1 mininum scan number is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->ms2MinScanNumEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS2 mininum scan number is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->pseudoMinPeakNumEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Pseudo mininum peak number is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->ms1EcscoreCutoffEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS1 ECScore cutoff is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->ms2EcscoreCutoffEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS2 ECScore cutoff is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->pseudoScoreCutoffEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Pseudo Score cutoff is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->ms1IntePccCutoffEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS1 seed envelope intensity cutoff is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->ms2IntePccCutoffEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MS2 seed envelope intensity cutoff is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->threadNumberEdit->text().toInt() > toppic::mem_check::getMaxThreads("topdia")) {
    int max_thread = toppic::mem_check::getMaxThreads("topdia");
    QMessageBox::StandardButton reply = QMessageBox::warning(this, tr("Warning"),
                         QString("Thread number is too large! Based on the memory size, up to %1 threads can run on this computer. Are you sure you want to proceed?").arg(max_thread).arg(max_thread),
                         QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::No) {
      return true;
    }
  }

  return false;
}

void TopDIADialog::updateMsg(std::string msg) {
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

void TopDIADialog::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait) {
    QCoreApplication::processEvents();
  }
}

