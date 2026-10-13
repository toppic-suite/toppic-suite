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

#include "gui/topfd/topfddialog.hpp"

#include <QCloseEvent>
#include <QDebug>
#include <QDesktopServices>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QFontDatabase>
#include <QMessageBox>
#include <QProcess>
#include <QRegularExpression>
#include <QRegularExpressionValidator>
#include <QScrollBar>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "common/util/file_util.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/version.hpp"
#include "gui/topfd/ui_topfddialog.h"
#include "gui/util/command.hpp"
#include "gui/util/gui_check.hpp"
#include "gui/util/gui_message.hpp"

TopFDDialog::TopFDDialog(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::TopFDDialog) {
  para_ptr_ = std::make_shared<toppic::TopfdPara>();
  ui->setupUi(this);
  std::string title = "TopFD v." + toppic::Version::getVersion();
  QString qstr = QString::fromStdString(title);
  this->setWindowTitle(qstr);
  lastDir_ = ".";
  ui->maxChargeEdit->setValidator(new QIntValidator(1, 100, this));
  ui->maxMassEdit->setValidator(new QIntValidator(1, 1000000, this));
  // 0.0001 to 0.1
  QRegularExpression rx1("0\\.0\\d{0,2}[1-9]|0\\.1");
  QRegularExpressionValidator* validator1 =
      new QRegularExpressionValidator(rx1, this);
  ui->mzErrorEdit->setValidator(validator1);
  // non-negative, the decimal point optional (the defaults are "3" and "1")
  QRegularExpression rx2("\\d{1,6}(\\.\\d{0,2})?");
  QRegularExpressionValidator* validator2 =
      new QRegularExpressionValidator(rx2, this);
  ui->ms1snRatioEdit->setValidator(validator2);
  ui->ms2snRatioEdit->setValidator(validator2);
  ui->splitRatioEdit->setValidator(
      new QDoubleValidator(0.01, 1000000, 2, this));
  ui->threadNumberEdit->setValidator(new QIntValidator(1, 1000, this));
  ui->minScanNumEdit->setValidator(new QIntValidator(1, 3, this));
  ui->ecscoreCutoffEdit->setValidator(new QDoubleValidator(0.0, 1.0, 4, this));
  ui->ms1EnvCnnCutoffEdit->setValidator(
      new QDoubleValidator(0.0, 1.0, 4, this));
  ui->ms2EnvCnnCutoffEdit->setValidator(
      new QDoubleValidator(0.0, 1.0, 4, this));
  ui->windowSizeEdit->setValidator(new QDoubleValidator(0.01, 10000, 2, this));
  QFont font;
  QFont outputFont;
#if defined(_WIN32) || defined(_WIN64) || defined(__MINGW32__) || \
    defined(__MINGW64__)
  font.setFamily(QStringLiteral("Calibri"));
  outputFont.setFamily(QStringLiteral("Consolas"));
#else
  const QString monoFamily =
      QFontDatabase::systemFont(QFontDatabase::FixedFont).family();
  font.setFamily(monoFamily);
  outputFont.setFamily(monoFamily);
#endif
  font.setPixelSize(12);
  outputFont.setPixelSize(12);
  QApplication::setFont(font);
  ui->outputTextBrowser->setFont(outputFont);
  TopFDDialog::on_defaultButton_clicked();
}

TopFDDialog::~TopFDDialog() {
  if (process_.state() != QProcess::NotRunning) {
    process_.kill();
  }
  delete ui;
}

void TopFDDialog::closeEvent(QCloseEvent* event) {
  if (process_.state() != QProcess::NotRunning) {
    if (!continueToClose()) {
      event->ignore();
      return;
    } else {
      process_.kill();
    }
  }
  event->accept();
  return;
}

void TopFDDialog::on_clearButton_clicked() {
  ui->listWidget->clear();
  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText(
      "Click the Start button to process the spectrum files.");
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
  ui->splitRatioEdit->setText(
      QString::number(para_ptr->getSplitIntensityRatio()));
  ui->windowSizeEdit->setText(QString::number(para_ptr->getPrecWindowWidth()));
  ui->threadNumberEdit->setText(QString::number(para_ptr->getThreadNum()));
  ui->ecscoreCutoffEdit->setText(
      QString::number(para_ptr->getMs1EcscoreCutoff()));
  ui->ms1EnvCnnCutoffEdit->setText(
      QString::number(para_ptr->getMs1EnvCnnScoreCutoff()));
  ui->ms2EnvCnnCutoffEdit->setText(
      QString::number(para_ptr->getMs2EnvCnnScoreCutoff()));
  ui->minScanNumEdit->setText(QString::number(para_ptr->getMs1MinScanNum()));
  ui->msDeconvCheckBox->setChecked(para_ptr->isSortUseMsDeconv());
  ui->missLevelOneCheckBox->setChecked(para_ptr->isMissingLevelOne());
  ui->fragNumFilteringCheckBox->setChecked(para_ptr->isAANumBasedFilter());
  ui->disableAdditionalFeatureSearchCheckBox->setChecked(
      !para_ptr->isSearchPrecWindow());
  ui->singleScanNoiseLevelCheckBox->setChecked(
      para_ptr->isUseSingleScanNoiseLevel());
  ui->noSqlCheckBox->setChecked(!para_ptr->isGeneSql());
  ui->ms1DeconvComboBox->setCurrentText(
      QString::fromStdString(para_ptr->getMs1DeconvMethod()));
  ui->ms2DeconvComboBox->setCurrentText(
      QString::fromStdString(para_ptr->getMs2DeconvMethod()));
  ui->textPeakListCheckBox->setChecked(para_ptr->isTextPeakList());

  ui->outputTextBrowser->clear();
  ui->outputTextBrowser->setText(
      "Click the Start button to process the spectrum files.");
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
  // a text peak list may have any extension (topfd --text-peak-list); an
  // MGF file of one spectrum is read as one
  QString filter = ui->textPeakListCheckBox->isChecked()
                       ? "Text peak lists (*.txt *.mgf);;All files (*)"
                       : "Spectra files (*.mzXML *.mzML *.mzxml *.mzml)";
  QStringList spfiles = QFileDialog::getOpenFileNames(
      this, "Select spectrum files", lastDir_, filter);
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
    // lastDir_ = s;
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
  QListWidgetItem* delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
}

void TopFDDialog::on_startButton_clicked() {
  if (checkError()) {
    return;
  }
  lockDialog();
  toppic::TopfdParaPtr paraPtr = this->getParaPtr();
  std::vector<std::string> specFileList = this->getSpecFileList();

  std::string cmd =
      toppic::command::geneTopfdCommand(para_ptr_, spec_file_lst_);
  QString q_cmd = QString::fromStdString(cmd);
  q_cmd = q_cmd.trimmed();
  QStringList cmd_list = q_cmd.split(" ");
  QString prog = cmd_list[0];
  cmd_list.removeFirst();

  // qDebug() << "start process ";
  process_.start(prog, cmd_list);
  process_.waitForStarted();
  // qDebug() << "start process finished";

  toppic::GuiMessage guiMsg;
  bool finish = false;
  while (!finish) {
    if (process_.state() == QProcess::NotRunning) {
      finish = true;
    }
    bool ready = process_.waitForReadyRead(100);
    if (ready || finish) {
      // qDebug() << "read finished";
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
        str = str + "ERROR Quit code: " + QString::number(process_.exitCode()) +
              ".\n";
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

void TopFDDialog::on_exitButton_clicked() { close(); }

bool TopFDDialog::continueToClose() {
  if (QMessageBox::question(
          this, tr("Quit"),
          tr("TopFD is still running. Are you sure you want to quit?"),
          QMessageBox::Yes | QMessageBox::No,
          QMessageBox::No) == QMessageBox::Yes) {
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
    ui->outputTextBrowser->setText(
        "Current directory " + QString::fromStdString(exe_dir) +
        " contains space and will cause errors in the program!");
  }
  para_ptr_->setResourceDir(toppic::file_util::getResourceDir(exe_dir));
  para_ptr_->setMaxCharge(std::stoi(ui->maxChargeEdit->text().toStdString()));
  para_ptr_->setMaxMass(std::stod(ui->maxMassEdit->text().toStdString()));
  para_ptr_->setMzError(std::stod(ui->mzErrorEdit->text().toStdString()));
  para_ptr_->setMsOneSnRatio(
      std::stod(ui->ms1snRatioEdit->text().toStdString()));
  para_ptr_->setMsTwoSnRatio(
      std::stod(ui->ms2snRatioEdit->text().toStdString()));
  para_ptr_->setSplitIntensityRatio(
      std::stod(ui->splitRatioEdit->text().toStdString()));
  para_ptr_->setPrecWindowWidth(
      std::stod(ui->windowSizeEdit->text().toStdString()));
  para_ptr_->setMissingLevelOne(ui->missLevelOneCheckBox->isChecked());
  para_ptr_->setMs1EcscoreCutoff(
      std::stod(ui->ecscoreCutoffEdit->text().toStdString()));
  para_ptr_->setMs1EnvCnnScoreCutoff(
      std::stod(ui->ms1EnvCnnCutoffEdit->text().toStdString()));
  para_ptr_->setMs2EnvCnnScoreCutoff(
      std::stod(ui->ms2EnvCnnCutoffEdit->text().toStdString()));
  para_ptr_->setMs1MinScanNum(
      std::stoi(ui->minScanNumEdit->text().toStdString()));
  para_ptr_->setThreadNum(
      std::stoi(ui->threadNumberEdit->text().toStdString()));
  para_ptr_->setDeconvScore(
      ui->msDeconvCheckBox->isChecked() ? "msdeconv" : "envcnn");
  para_ptr_->setActivation(ui->activationComboBox->currentText().toStdString());
  para_ptr_->setAANumBasedFilter(ui->fragNumFilteringCheckBox->isChecked());

  para_ptr_->setSearchPrecWindow(
      (ui->disableAdditionalFeatureSearchCheckBox->isChecked()));
  para_ptr_->setUseSingleScanNoiseLevel(
      ui->singleScanNoiseLevelCheckBox->isChecked());
  para_ptr_->setGeneSql(!ui->noSqlCheckBox->isChecked());
  para_ptr_->setMs1DeconvMethod(
      ui->ms1DeconvComboBox->currentText().toStdString());
  para_ptr_->setMs2DeconvMethod(
      ui->ms2DeconvComboBox->currentText().toStdString());
  para_ptr_->setTextPeakList(ui->textPeakListCheckBox->isChecked());

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
  ui->splitRatioEdit->setEnabled(false);
  ui->threadNumberEdit->setEnabled(false);
  ui->minScanNumEdit->setEnabled(false);
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->missLevelOneCheckBox->setEnabled(false);
  ui->windowSizeEdit->setEnabled(false);
  ui->outputButton->setEnabled(false);
  ui->msDeconvCheckBox->setEnabled(false);
  ui->activationComboBox->setEnabled(false);
  ui->fragNumFilteringCheckBox->setEnabled(false);
  ui->ecscoreCutoffEdit->setEnabled(false);
  ui->ms1EnvCnnCutoffEdit->setEnabled(false);
  ui->ms2EnvCnnCutoffEdit->setEnabled(false);
  ui->disableAdditionalFeatureSearchCheckBox->setEnabled(false);
  ui->singleScanNoiseLevelCheckBox->setEnabled(false);
  ui->noSqlCheckBox->setEnabled(false);
  ui->ms1DeconvComboBox->setEnabled(false);
  ui->ms2DeconvComboBox->setEnabled(false);
  ui->textPeakListCheckBox->setEnabled(false);
}

void TopFDDialog::unlockDialog() {
  ui->addButton->setEnabled(true);
  ui->delButton->setEnabled(true);
  ui->maxChargeEdit->setEnabled(true);
  ui->maxMassEdit->setEnabled(true);
  ui->mzErrorEdit->setEnabled(true);
  ui->ms1snRatioEdit->setEnabled(true);
  ui->ms2snRatioEdit->setEnabled(true);
  ui->splitRatioEdit->setEnabled(true);
  ui->threadNumberEdit->setEnabled(true);
  ui->minScanNumEdit->setEnabled(true);
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->missLevelOneCheckBox->setEnabled(true);
  ui->windowSizeEdit->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);
  ui->msDeconvCheckBox->setEnabled(true);
  ui->activationComboBox->setEnabled(true);
  ui->fragNumFilteringCheckBox->setEnabled(true);
  ui->ecscoreCutoffEdit->setEnabled(true);
  ui->ms1EnvCnnCutoffEdit->setEnabled(true);
  ui->ms2EnvCnnCutoffEdit->setEnabled(true);
  ui->disableAdditionalFeatureSearchCheckBox->setEnabled(true);
  ui->singleScanNoiseLevelCheckBox->setEnabled(true);
  ui->noSqlCheckBox->setEnabled(true);
  ui->ms1DeconvComboBox->setEnabled(true);
  ui->ms2DeconvComboBox->setEnabled(true);
  ui->textPeakListCheckBox->setEnabled(true);
}

bool TopFDDialog::checkError() {
  if (ui->listWidget->count() == 0) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Please add at least one spectrum file!"),
                         QMessageBox::Yes);
    return true;
  }
  for (int i = 0; i < ui->listWidget->count(); i++) {
    if (toppic::gui_check::pathWithSpace(
            this, ui->listWidget->item(i)->text(), tr("spectrum file"))) {
      return true;
    }
  }

  const std::vector<std::pair<QLineEdit*, QString>> edits = {
      {ui->maxChargeEdit, tr("Maximum charge")},
      {ui->maxMassEdit, tr("Maximum mass")},
      {ui->mzErrorEdit, tr("M/z error")},
      {ui->windowSizeEdit, tr("Default precursor window")},
      {ui->ms1snRatioEdit, tr("MS1 S/N ratio")},
      {ui->ms2snRatioEdit, tr("MS2 S/N ratio")},
      {ui->threadNumberEdit, tr("Thread number")},
      {ui->ms1EnvCnnCutoffEdit, tr("MS1 EnvCNN cutoff")},
      {ui->ms2EnvCnnCutoffEdit, tr("MS2 EnvCNN cutoff")},
      {ui->ecscoreCutoffEdit, tr("ECScore cutoff")},
      {ui->splitRatioEdit, tr("Feature split ratio")},
      {ui->minScanNumEdit, tr("Feature minimum scan number")}};
  for (const auto& [edit, name] : edits) {
    if (toppic::gui_check::invalidEdit(this, edit, name)) {
      return true;
    }
  }

  if (ui->textPeakListCheckBox->isChecked()) {
    // topfd --text-peak-list reads every input file as a peak list; an
    // mzML/mzXML file has no line that starts with a number and would give
    // an empty spectrum
    for (int i = 0; i < ui->listWidget->count(); i++) {
      QString spfile = ui->listWidget->item(i)->text();
      if (spfile.endsWith(".mzML", Qt::CaseInsensitive) ||
          spfile.endsWith(".mzXML", Qt::CaseInsensitive)) {
        QMessageBox::warning(
            this, tr("Warning"),
            tr("Text peak list input is checked, but %1 is an mzML or mzXML "
               "file. Remove it or uncheck the option.")
                .arg(spfile),
            QMessageBox::Yes);
        return true;
      }
    }
    // the FLASHDeconv process reads mzML files only; topfd refuses the pair
    if (ui->ms2DeconvComboBox->currentText() == "flashdeconv") {
      QMessageBox::warning(this, tr("Warning"),
                           tr("The MS2 deconvolution method flashdeconv "
                              "cannot be used with text peak list input!"),
                           QMessageBox::Yes);
      return true;
    }
  }

  if (ui->threadNumberEdit->text().toInt() >
      toppic::mem_check::getMaxThreads("topfd")) {
    int max_thread = toppic::mem_check::getMaxThreads("topfd");
    QMessageBox::StandardButton reply = QMessageBox::warning(
        this, tr("Warning"),
        QString("Thread number is too large! Based on the memory size, up to "
                "%1 threads can run on this computer. Are you sure you want to "
                "proceed?")
            .arg(max_thread),
        QMessageBox::Yes | QMessageBox::No);
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
