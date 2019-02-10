//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QCloseEvent>
#include <QCoreApplication>
#include <QToolTip>
#include <QDesktopServices>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "toppicwindow.h"
#include "ui_toppicwindow.h"
#include "threadtoppic.h"

toppicWindow::toppicWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::toppicWindow) {
      ui->setupUi(this);
      lastDir_ = ".";
      QRegExp rx1("^\\d{1,8}\\.\\d{0,2}$");
      QRegExpValidator *validator1 = new QRegExpValidator(rx1, this);
      ui->maxModEdit->setValidator(validator1);
      ui->cutoffSpectralValueEdit->setValidator(validator1);
      ui->cutoffProteoformValueEdit->setValidator(validator1);
      ui->numCombinedEdit->setValidator(new QIntValidator(0, 2147483647, this));
      QRegExp rx2("^0\\.\\d{0,2}|1.00$");
      QRegExpValidator *validator2 = new QRegExpValidator(rx2, this);
      ui->miscoreThresholdEdit->setValidator(validator2);
      ui->threadNumberEdit->setValidator(new QIntValidator(0, 2147483647, this));
      ui->errorToleranceEdit->setValidator(new QIntValidator(0, 2147483647, this));
      QRegExp rx3("^-?\\d{1,8}\\.\\d{0,2}$");
      QRegExpValidator *validator3 = new QRegExpValidator(rx3, this);
      ui->minModEdit->setValidator(validator3);
      QFont font;
      QFont fontTable;
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
      font.setFamily(QStringLiteral("Courier New"));
      fontTable.setFamily(QStringLiteral("Courier New"));
#else
      font.setFamily(QStringLiteral("Monospace"));
      fontTable.setFamily(QStringLiteral("Monospace"));
#endif
      font.setPointSize(12);
      QApplication::setFont(font);
      ui->outputTextBrowser->setFont(font);
      thread_ = new threadtoppic(this);
      showInfo = "";
      setToolTip("");
      setToolTipDuration(100);

      fontTable.setPointSize(9);
      ui->listWidget->setFont(fontTable);

      on_clearButton_clicked();
      on_defaultButton_clicked();
      ui->tabWidget->setCurrentIndex(0);
    }

toppicWindow::~toppicWindow() {
  delete ui;
}

void toppicWindow::initArguments() {
  arguments_["oriDatabaseFileName"]="";
  arguments_["databaseFileName"] = "";
  arguments_["databaseBlockSize"] = "1000000";
  arguments_["combinedOutputName"] = "";
  arguments_["activation"] = "FILE";
  arguments_["searchType"] = "TARGET";
  arguments_["fixedMod"] = "";
  arguments_["ptmNumber"] = "1";
  arguments_["errorTolerance"] = "15";
  arguments_["cutoffSpectralType"] = "EVALUE";
  arguments_["cutoffSpectralValue"] = "0.01";
  arguments_["cutoffProteoformType"] = "EVALUE";
  arguments_["cutoffProteoformValue"] = "0.01";
  arguments_["allowProtMod"] = "NONE,NME,NME_ACETYLATION,M_ACETYLATION";
  arguments_["numOfTopPrsms"] = "1";
  arguments_["maxPtmMass"] = "500";
  arguments_["minPtmMass"] = "-500";
  arguments_["useLookupTable"] = "false";
  arguments_["executiveDir"] = ".";
  arguments_["resourceDir"] = "";
  arguments_["keepTempFiles"] = "false";
  arguments_["localThreshold"] = "0.45";
  arguments_["groupSpectrumNumber"] = "1";
  arguments_["filteringResultNumber"] = "20";
  arguments_["residueModFileName"] = "";
  arguments_["threadNumber"] = "1";
  arguments_["useFeatureFile"] = "true";
  arguments_["skipList"] = "";
}

void toppicWindow::on_clearButton_clicked() {
  ui->databaseFileEdit->clear();
  ui->listWidget->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
  ui->combinedOutputEdit->setText("");
  ui->combinedOutputEdit->setEnabled(false);
  ui->outputButton->setEnabled(false);
}

void toppicWindow::on_defaultButton_clicked() {
  ui->combinedOutputEdit->setText("");
  ui->fixedModFileEdit->clear();
  ui->errorToleranceEdit->setText("15");
  ui->maxModEdit->setText("500");
  ui->minModEdit->setText("-500");
  ui->cutoffSpectralValueEdit->setText("0.01");
  ui->cutoffProteoformValueEdit->setText("0.01");
  ui->numCombinedEdit->setText("1");
  ui->miscoreThresholdEdit->setText("0.45");
  ui->threadNumberEdit->setText("1");
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
  ui->fixedModComboBox->setCurrentIndex(0);
  on_fixedModComboBox_currentIndexChanged(0);
  ui->activationComboBox->setCurrentIndex(0);
  ui->cutoffSpectralTypeComboBox->setCurrentIndex(0);
  ui->cutoffProteoformTypeComboBox->setCurrentIndex(0);
  ui->numModComboBox->setCurrentIndex(1);
  on_numModComboBox_currentIndexChanged(1);
  ui->NONECheckBox->setChecked(true);
  ui->NMECheckBox->setChecked(true);
  ui->NMEACCheckBox->setChecked(true);
  ui->MACCheckBox->setChecked(true);
  ui->decoyCheckBox->setChecked(false);
  ui->lookupTableCheckBox->setChecked(false);
  ui->topfdFeatureCheckBox->setChecked(false);
}

void toppicWindow::updatedir(QString s) {
  if (!s.isEmpty()) {
    lastDir_ = s;
  }
}

// single file
void toppicWindow::toppicWindow::on_databaseFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a protein database file",
      lastDir_,
      "Database files(*.fasta *.fa)");
  updatedir(s);
  ui->databaseFileEdit->setText(s);
}

void toppicWindow::on_fixedModFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a fixed modification file",
      lastDir_,
      "Modification files(*.txt);;All files(*.*)");
  updatedir(s);
  ui->fixedModFileEdit->setText(s);
}

void toppicWindow::on_modFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a modification file for PTM localization",
      lastDir_,
      "Modification files(*.txt);;All files(*.*)");
  updatedir(s);
  ui->modFileEdit->setText(s);
}

void toppicWindow::on_startButton_clicked() {
  
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

void toppicWindow::on_outputButton_clicked() {
  std::vector<std::string> spec_file_lst = this->getSpecFileList();
  if (spec_file_lst.size() > 0) {
    fs::path full_path(spec_file_lst[0].c_str());
    QString outPath = full_path.remove_filename().string().c_str();
    QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
  }
}

std::map<std::string, std::string> toppicWindow::getArguments() {
  QString path = QCoreApplication::applicationFilePath();
  std::string exe_dir = prot::file_util::getExecutiveDir(path.toStdString());
  arguments_["executiveDir"] = exe_dir;
  arguments_["resourceDir"] = arguments_["executiveDir"] + prot::file_util::getFileSeparator() + prot::file_util::getResourceDirName();
  arguments_["oriDatabaseFileName"] = ui->databaseFileEdit->text().toStdString();
  
  arguments_["combinedOutputName"] = ui->combinedOutputEdit->text().trimmed().toStdString();
  arguments_["databaseBlockSize"] = "1000000";
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
  arguments_["ptmNumber"] = ui->numModComboBox->currentText().toStdString();
  arguments_["errorTolerance"] = ui->errorToleranceEdit->text().toStdString();
  arguments_["cutoffSpectralType"] = ui->cutoffSpectralTypeComboBox->currentText().toStdString();
  arguments_["cutoffSpectralValue"] = ui->cutoffSpectralValueEdit->text().toStdString();
  arguments_["cutoffProteoformType"] = ui->cutoffProteoformTypeComboBox->currentText().toStdString();
  arguments_["cutoffProteoformValue"] = ui->cutoffProteoformValueEdit->text().toStdString();
  arguments_["allowProtMod"] = "";
  if (ui->NONECheckBox->isChecked()) {
    arguments_["allowProtMod"] = arguments_["allowProtMod"] + ",NONE";
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
  if (arguments_["allowProtMod"] != "") {
    arguments_["allowProtMod"] = arguments_["allowProtMod"].substr(1);
  }
  arguments_["numOfTopPrsms"] = "1";
  arguments_["maxPtmMass"] = ui->maxModEdit->text().toStdString();
  arguments_["minPtmMass"] = ui->minModEdit->text().toStdString();
  if (ui->lookupTableCheckBox->isChecked()) {
    arguments_["useLookupTable"] = "true";
  } else {
    arguments_["useLookupTable"] = "false";
  }
  arguments_["keepTempFiles"] = "false";   // default
  arguments_["localThreshold"] = ui->miscoreThresholdEdit->text().toStdString();
  arguments_["groupSpectrumNumber"] = ui->numCombinedEdit->text().toStdString();
  arguments_["filteringResultNumber"] = "20";  // default
  arguments_["residueModFileName"] = ui->modFileEdit->text().toStdString();
  arguments_["threadNumber"] = ui->threadNumberEdit->text().toStdString();
  if (ui->topfdFeatureCheckBox->isChecked()) {
    arguments_["useFeatureFile"] = "false";
  } else {
    arguments_["useFeatureFile"] = "true";
  }
  //showArguments();
  return arguments_;
}

std::vector<std::string> toppicWindow::getSpecFileList() {
  spec_file_lst_.clear();
  for (int i = 0; i < ui->listWidget->count(); i++) {
    spec_file_lst_.push_back(ui->listWidget->item(i)->text().toStdString());
  }
  return spec_file_lst_;
}

void toppicWindow::on_addButton_clicked() {
  QStringList spfiles = QFileDialog::getOpenFileNames(
      this,
      "Select deconvoluted spectrum files",
      lastDir_,
      "Spectrum files(*ms2.msalign)");

  for (int i = 0; i < spfiles.size(); i++) {
    QString spfile = spfiles.at(i);
    updatedir(spfile);
    if (ableToAdd(spfile)) {
      ui->listWidget->addItem(new QListWidgetItem(spfile));
    }
  }
  
  if (ui->listWidget->count() > 1) {
    ui->combinedOutputEdit->setEnabled(true);
  }
}

bool toppicWindow::ableToAdd(QString spfile) {
  bool able = true;
  if (spfile != "") {
    if (spfile.toStdString().length() > 200) {
      QMessageBox::warning(this, tr("Warning"),
                           tr("The deconvoluted spectrum file path is too long!"),
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

void toppicWindow::on_delButton_clicked() {
  QListWidgetItem *delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
  if (ui->listWidget->count() < 2) {
    ui->combinedOutputEdit->setEnabled(false);
  }
}

void toppicWindow::lockDialog() {
  ui->databaseFileButton->setEnabled(false);
  ui->modFileButton->setEnabled(false);
  ui->databaseFileEdit->setEnabled(false);
  ui->combinedOutputEdit->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  ui->fixedModFileEdit->setEnabled(false);
  ui->fixedModFileButton->setEnabled(false);
  ui->errorToleranceEdit->setEnabled(false);
  ui->maxModEdit->setEnabled(false);
  ui->minModEdit->setEnabled(false);
  ui->cutoffSpectralValueEdit->setEnabled(false);
  ui->cutoffProteoformValueEdit->setEnabled(false);
  ui->numCombinedEdit->setEnabled(false);
  ui->modFileEdit->setEnabled(false);
  ui->miscoreThresholdEdit->setEnabled(false);
  ui->threadNumberEdit->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  ui->activationComboBox->setEnabled(false);
  ui->cutoffSpectralTypeComboBox->setEnabled(false);
  ui->cutoffProteoformTypeComboBox->setEnabled(false);
  ui->numModComboBox->setEnabled(false);
  ui->NONECheckBox->setEnabled(false);
  ui->NMECheckBox->setEnabled(false);
  ui->NMEACCheckBox->setEnabled(false);
  ui->MACCheckBox->setEnabled(false);
  ui->decoyCheckBox->setEnabled(false);
  ui->lookupTableCheckBox->setEnabled(false);
  ui->topfdFeatureCheckBox->setEnabled(false);
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->outputButton->setEnabled(false);
  ui->addButton->setEnabled(false);
  ui->delButton->setEnabled(false);
}

void toppicWindow::unlockDialog() {
  ui->databaseFileButton->setEnabled(true);
  ui->modFileButton->setEnabled(true);
  ui->databaseFileEdit->setEnabled(true);
  if (ui->listWidget->count() > 1) {
    ui->combinedOutputEdit->setEnabled(true);
  }
  ui->fixedModFileEdit->setEnabled(true);
  ui->errorToleranceEdit->setEnabled(true);
  ui->maxModEdit->setEnabled(true);
  ui->minModEdit->setEnabled(true);
  ui->cutoffSpectralValueEdit->setEnabled(true);
  ui->cutoffProteoformValueEdit->setEnabled(true);
  ui->numCombinedEdit->setEnabled(true);
  ui->modFileEdit->setEnabled(true);
  ui->miscoreThresholdEdit->setEnabled(true);
  ui->threadNumberEdit->setEnabled(true);
  ui->fixedModComboBox->setEnabled(true);
  on_fixedModComboBox_currentIndexChanged(ui->fixedModComboBox->currentIndex());
  ui->activationComboBox->setEnabled(true);
  ui->cutoffSpectralTypeComboBox->setEnabled(true);
  ui->cutoffProteoformTypeComboBox->setEnabled(true);
  ui->numModComboBox->setEnabled(true);
  ui->NONECheckBox->setEnabled(true);
  ui->NMECheckBox->setEnabled(true);
  ui->NMEACCheckBox->setEnabled(true);
  ui->MACCheckBox->setEnabled(true);
  ui->decoyCheckBox->setEnabled(true);
  ui->lookupTableCheckBox->setEnabled(true);
  ui->topfdFeatureCheckBox->setEnabled(true);
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);
  ui->addButton->setEnabled(true);
  ui->delButton->setEnabled(true);
}

bool toppicWindow::checkError() {
  if (ui->databaseFileEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Please select a protein database file!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->databaseFileEdit->text().toStdString().length() > 200) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("The protein database file path is too long!"),
                         QMessageBox::Yes);
    return true;
  }

  for (int i = 0; i < ui->listWidget->count(); i++) {
    if (ui->listWidget->item(i)->text().toStdString().length() > 200) {
      QMessageBox::warning(this, tr("Warning"),
                           tr("The sepctrum file path is too long!"),
                           QMessageBox::Yes);
      return true;
    }
  }

  if (ui->combinedOutputEdit->text().toStdString().length() > 200) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("The combined output filename is too long!"),
                         QMessageBox::Yes);
    return true;
  }

  QString currentText = ui->errorToleranceEdit->text();
  if (ui->lookupTableCheckBox->isChecked() && currentText != "5" && currentText != "10" && currentText != "15") {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an error tolerance other than 5, 10, and 15 ppm, the checkbox \"Lookup table for E-value computation\" should be not selected!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->fixedModFileEdit->text().isEmpty() && ui->fixedModComboBox->currentIndex() == 3) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Please select a fixed modification file!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->errorToleranceEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Error tolerance is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->maxModEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Maximum mass shift is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->minModEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Minimum mass shift is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->cutoffSpectralValueEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Spectrum-level cutoff value is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->cutoffProteoformValueEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Proteoform-level cutoff value is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->numCombinedEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Number of combined spectra is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->miscoreThresholdEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("MIScore threshold is empty!"),
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

void toppicWindow::updateMsg(std::string msg) {
  showInfo = msg.c_str();
  ui->outputTextBrowser->setText(showInfo);
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
  QString info = msg.c_str();
}

void toppicWindow::showArguments() {
  QMessageBox::warning(0, "Arguments", ("executiveDir:" + arguments_["executiveDir"] +
                                        "\nresourceDir:" + arguments_["resourceDir"] +
                                        "\noriDatabaseFileName:" + arguments_["oriDatabaseFileName"] +
                                        "\ndatabaseFileName:" + arguments_["databaseFileName"] +
                                        "\ndatabaseBlockSize:" + arguments_["databaseBlockSize"] +
                                        "\ncombinedOutputName:" + arguments_["combinedOutputName"] +
                                        "\nactivation:" + arguments_["activation"] +
                                        "\nsearchType:" + arguments_["searchType"] +
                                        "\nfixedMod:" + arguments_["fixedMod"] +
                                        "\nptmNumber:" + arguments_["ptmNumber"] +
                                        "\nerrorTolerance:" + arguments_["errorTolerance"] +
                                        "\ncutoffSpectralType:" + arguments_["cutoffSpectralType"] +
                                        "\ncutoffSpectralValue:" + arguments_["cutoffSpectralValue"] +
                                        "\ncutoffProteoformType:" + arguments_["cutoffProteoformType"] +
                                        "\ncutoffProteoformValue:" + arguments_["cutoffProteoformValue"] +
                                        "\nallowProtMod:" + arguments_["allowProtMod"] +
                                        "\nnumOfTopPrsms:" + arguments_["numOfTopPrsms"] +
                                        "\nmaxPtmMass:" + arguments_["maxPtmMass"] +
                                        "\nminPtmMass:" + arguments_["minPtmMass"] +
                                        "\nuseLookupTable:" + arguments_["useLookupTable"] +
                                        "\nkeepTempFiles:" + arguments_["keepTempFiles"] +
                                        "\nlocalThreshold:" + arguments_["localThreshold"] +
                                        "\ngroupSpectrumNumber:" + arguments_["groupSpectrumNumber"] +
                                        "\nfilteringResultNumber:" + arguments_["filteringResultNumber"] +
                                        "\nresidueModFileName:" + arguments_["residueModFileName"] +
                                        "\nthreadNumber:" + arguments_["threadNumber"] +
                                        "\nuseFeatureFile:" + arguments_["useFeatureFile"] +
                                        "\nskipList:" + arguments_["skipList"]).c_str(), QMessageBox::Yes);
}

void toppicWindow::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait)
    QCoreApplication::processEvents();
}

void toppicWindow::on_exitButton_clicked() {
  close();
}

void toppicWindow::on_fixedModComboBox_currentIndexChanged(int index) {
  if (index == 3) {
    ui->fixedModFileEdit->setEnabled(true);
    ui->fixedModFileButton->setEnabled(true);
  } else {
    ui->fixedModFileEdit->setEnabled(false);
    ui->fixedModFileButton->setEnabled(false);
  }
}

void toppicWindow::on_numModComboBox_currentIndexChanged(int index) {
  if (index == 0) {
    ui->modFileEdit->setEnabled(false);
    ui->modFileButton->setEnabled(false);
    ui->maxModEdit->setEnabled(false);
    ui->minModEdit->setEnabled(false);
    ui->miscoreThresholdEdit->setEnabled(false);
  } else {
    ui->modFileEdit->setEnabled(true);
    ui->modFileButton->setEnabled(true);
    ui->maxModEdit->setEnabled(true);
    ui->minModEdit->setEnabled(true);
    ui->miscoreThresholdEdit->setEnabled(true);
  }
}

void toppicWindow::on_errorToleranceEdit_textChanged(QString string) {
  QString currentText = ui->errorToleranceEdit->text();
  if (ui->lookupTableCheckBox->isChecked() && currentText != "5" && currentText != "10" && currentText != "15" && currentText != "1") {
    QMessageBox::warning(this, tr("Warning"),
                         tr("When the checkbox \"Lookup table for E-value computation\" is checked, only three error tolerance values 5, 10, and 15 ppm can be used!"),
                         QMessageBox::Yes);
    ui->errorToleranceEdit->setText("15");
  }
}

void toppicWindow::on_lookupTableCheckBox_clicked(bool checked) {
  QString currentText = ui->errorToleranceEdit->text();
  if (checked && currentText != "5" && currentText != "10" && currentText != "15") {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an error tolerance other than 5, 10, and 15 ppm, the checkbox \"Lookup table for E-value computation\" should not be checked!"),
                         QMessageBox::Yes);
    ui->lookupTableCheckBox->setChecked(true);
  }
}

bool toppicWindow::nterminalerror() {
  if (ui->NONECheckBox->isChecked() || ui->NMECheckBox->isChecked() || ui->NMEACCheckBox->isChecked() || ui->MACCheckBox->isChecked()) {
    return false;
  } else {
    QMessageBox::warning(this, tr("Warning"),
                         tr("At least one N-terminal form should be selected!"),
                         QMessageBox::Yes);
    return true;
  }
}

void toppicWindow::on_NONECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NONECheckBox->setChecked(true);
  }
}

void toppicWindow::on_NMECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMECheckBox->setChecked(true);
  }
}

void toppicWindow::on_NMEACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMEACCheckBox->setChecked(true);
  }
}

void toppicWindow::on_MACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->MACCheckBox->setChecked(true);
  }
}

void toppicWindow::on_cutoffSpectralTypeComboBox_currentIndexChanged(int index) {
  if (index == 1 && !ui->decoyCheckBox->isChecked()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an FDR cutoff, the checkbox \"decoy database\" should be checked!"),
                         QMessageBox::Yes);
    ui->cutoffSpectralTypeComboBox->setCurrentIndex(0);
  }
}

void toppicWindow::on_cutoffProteoformTypeComboBox_currentIndexChanged(int index) {
  if (index == 1 && !ui->decoyCheckBox->isChecked()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an FDR cutoff, the checkbox \"decoy database\" should be checked!"),
                         QMessageBox::Yes);
    ui->cutoffProteoformTypeComboBox->setCurrentIndex(0);
  }
}

void toppicWindow::on_decoyCheckBox_clicked(bool checked) {
  if (!checked && (ui->cutoffSpectralTypeComboBox->currentIndex() > 0 || ui->cutoffProteoformTypeComboBox->currentIndex() > 0)) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Because an FDR cutoff is selected, the checkbox \"decoy database\" cannot be unchecked."),
                         QMessageBox::Yes);
    ui->decoyCheckBox->setChecked(true);
  }
}

void toppicWindow::closeEvent(QCloseEvent *event) {
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
bool toppicWindow::continueToClose() {
  if (QMessageBox::question(this,
                            tr("Quit"),
                            tr("TopPIC is still running. Are you sure you want to quit?"),
                            QMessageBox::Yes | QMessageBox::No,
                            QMessageBox::No)
      == QMessageBox::Yes) {
    return true;
  } else {
    return false;
  }
}

bool toppicWindow::event(QEvent *event) {
  if (event->type() == QEvent::ToolTip) {
    QHelpEvent *helpEvent = static_cast<QHelpEvent *>(event);
    if (QRect(800, 230, 60, 60).contains(helpEvent->pos()) && ui->tabWidget->currentIndex() == 1) {
      QToolTip::showText(helpEvent->globalPos(), "To use an error tolerance other than \n5, 10, and 15 ppm, the checkbox \n\"Lookup table for E-value computation\" should not be selected!");
    } else {
      QToolTip::hideText();
      event->ignore();
    }
    return true;
  }
  return QWidget::event(event);
}
