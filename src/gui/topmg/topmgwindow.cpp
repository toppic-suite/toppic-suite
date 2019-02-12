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

#include "topmgwindow.h"
#include "ui_topmgwindow.h"
#include "threadtopmg.h"

topmgWindow::topmgWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::topmgWindow) {
      ui->setupUi(this);
      lastDir_ = ".";
      QRegExp rx1("^\\d{1,8}\\.\\d{0,2}$");
      QRegExpValidator *validator1 = new QRegExpValidator(rx1, this);
      ui->maxModEdit->setValidator(validator1);
      ui->cutoffSpectralValueEdit->setValidator(validator1);
      ui->cutoffProteoformValueEdit->setValidator(validator1);
      ui->threadNumberEdit->setValidator(new QIntValidator(0, 2147483647, this));
      ui->errorToleranceEdit->setValidator(new QIntValidator(0, 2147483647, this));
      QFont font;
      QFont fontTable;
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
      font.setFamily(QStringLiteral("Calibri"));
      fontTable.setFamily(QStringLiteral("Calibri"));
#else
      font.setFamily(QStringLiteral("Monospace"));
      fontTable.setFamily(QStringLiteral("Monospace"));
#endif
      font.setPixelSize(12);
      QApplication::setFont(font);
      ui->outputTextBrowser->setFont(font);
      thread_ = new threadtopmg(this);
      showInfo = "";
      setToolTip("");
      setToolTipDuration(100);

      fontTable.setPixelSize(9);
      ui->listWidget->setFont(fontTable);

      on_clearButton_clicked();
      on_defaultButton_clicked();
      ui->tabWidget->setCurrentIndex(0);
    }

topmgWindow::~topmgWindow() {
  delete ui;
}

void topmgWindow::initArguments() {
  arguments_["oriDatabaseFileName"]="";
  arguments_["databaseFileName"] = "";
  arguments_["databaseBlockSize"] = "1000000";
  arguments_["spectrumFileName"] = "";
  arguments_["combinedOutputName"] = "";
  arguments_["activation"] = "FILE";
  arguments_["searchType"] = "TARGET";
  arguments_["fixedMod"] = "";
  arguments_["ptmNumber"] = "0";
  arguments_["errorTolerance"] = "15";
  arguments_["cutoffSpectralType"] = "EVALUE";
  arguments_["cutoffSpectralValue"] = "0.01";
  arguments_["cutoffProteoformType"] = "EVALUE";
  arguments_["cutoffProteoformValue"] = "0.01";
  arguments_["allowProtMod"] = "NONE,NME,NME_ACETYLATION,M_ACETYLATION";
  arguments_["numOfTopPrsms"] = "1";
  arguments_["maxPtmMass"] = "500";
  arguments_["executiveDir"] = ".";
  arguments_["resourceDir"] = "";
  arguments_["keepTempFiles"] = "false";
  arguments_["groupSpectrumNumber"] = "1";
  arguments_["filteringResultNumber"] = "20";
  arguments_["varModFileName"] = "";
  arguments_["threadNumber"] = "1";
  arguments_["useFeatureFile"] = "true";
  arguments_["skipList"] = "";
  arguments_["proteoGraphGap"] = "40";
  arguments_["useAsfDiag"] = "false";
  arguments_["varPtmNumber"] = "10";
  arguments_["varPtmNumInGap"] = "5";
}

void topmgWindow::on_clearButton_clicked() {
  ui->databaseFileEdit->clear();
  ui->listWidget->clear();
  ui->combinedOutputEdit->clear();
  ui->combinedOutputEdit->setEnabled(false);
  ui->skipListEdit->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
}

void topmgWindow::on_defaultButton_clicked() {
  ui->fixedModFileEdit->clear();
  ui->errorToleranceEdit->setText("15");
  ui->maxModEdit->setText("500");
  ui->cutoffSpectralValueEdit->setText("0.01");
  ui->cutoffProteoformValueEdit->setText("0.01");
  ui->threadNumberEdit->setText("1");
  ui->fixedModComboBox->setCurrentIndex(0);
  ui->fixedModFileEdit->setEnabled(false);
  ui->fixedModFileButton->setEnabled(false);
  ui->activationComboBox->setCurrentIndex(0);
  ui->cutoffSpectralTypeComboBox->setCurrentIndex(0);
  ui->cutoffProteoformTypeComboBox->setCurrentIndex(0);
  ui->numModComboBox->setCurrentIndex(4);
  ui->numUnknownShiftComboBox->setCurrentIndex(0);
  ui->NONECheckBox->setChecked(true);
  ui->NMECheckBox->setChecked(true);
  ui->NMEACCheckBox->setChecked(true);
  ui->MACCheckBox->setChecked(true);
  ui->decoyCheckBox->setChecked(false);
  ui->topfdFeatureCheckBox->setChecked(false);
  ui->asfDiagCheckBox->setChecked(false);
}

void topmgWindow::updatedir(QString s) {
  if (!s.isEmpty()) {
    lastDir_ = s;
  }
}

void topmgWindow::topmgWindow::on_databaseFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a protein database file",
      lastDir_,
      "Database files(*.fasta *.fa)");
  updatedir(s);
  ui->databaseFileEdit->setText(s);
}

void topmgWindow::on_fixedModFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a fixed modification file",
      lastDir_,
      "Modification files(*.txt);;All files(*.*)");
  updatedir(s);
  ui->fixedModFileEdit->setText(s);
}

void topmgWindow::on_modFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a modification file for variable PTMs",
      lastDir_,
      "Modification files(*.txt);;All files(*.*)");
  updatedir(s);
  ui->modFileEdit->setText(s);
}

void topmgWindow::on_skipListButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a modification file for skip list",
      lastDir_,
      "Text files(*.txt);;All files(*.*)");
  updatedir(s);
  ui->skipListEdit->setText(s);
}

void topmgWindow::on_startButton_clicked() {
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

void topmgWindow::on_outputButton_clicked() {
  std::vector<std::string> spec_file_lst = this->getSpecFileList();
  if (spec_file_lst.size() > 0) {
    fs::path full_path(spec_file_lst[0].c_str());
    QString outPath = full_path.remove_filename().string().c_str();
    QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
  }
}

std::map<std::string, std::string> topmgWindow::getArguments() {
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
  arguments_["varPtmNumber"] = ui->numModComboBox->currentText().toStdString();
  arguments_["ptmNumber"] = ui->numUnknownShiftComboBox->currentText().toStdString();
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
  arguments_["keepTempFiles"] = "false";   // default
  arguments_["filteringResultNumber"] = "20";  // default
  arguments_["useGf"] = "false";  // default
  arguments_["groupSpectrumNumber"] = "1";  // default
  arguments_["skipList"] = ui->skipListEdit->text().toStdString();
  arguments_["proteoGraphGap"] = "40";  // default
  arguments_["useAsfDiag"] = "false";  // default
  arguments_["varPtmNumber"] = ui->numModComboBox->currentText().toStdString();
  arguments_["ptmNumber"] = ui->numUnknownShiftComboBox->currentText().toStdString();
  arguments_["varPtmNumInGap"] = "5";  // default
  arguments_["varModFileName"] = ui->modFileEdit->text().toStdString();
  arguments_["threadNumber"] = ui->threadNumberEdit->text().toStdString();
  if (ui->topfdFeatureCheckBox->isChecked()) {
    arguments_["useFeatureFile"] = "false";
  } else {
    arguments_["useFeatureFile"] = "true";
  }
  if (ui->asfDiagCheckBox->isChecked()) {
    arguments_["useAsfDiag"] = "true";
  } else {
    arguments_["useAsfDiag"] = "false";
  }
  //showArguments();
  return arguments_;
}

std::vector<std::string> topmgWindow::getSpecFileList() {
  spec_file_lst_.clear();
  for (int i = 0; i < ui->listWidget->count(); i++) {
    spec_file_lst_.push_back(ui->listWidget->item(i)->text().toStdString());
  }
  return spec_file_lst_;
}

void topmgWindow::on_addButton_clicked() {
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

bool topmgWindow::ableToAdd(QString spfile) {
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

void topmgWindow::on_delButton_clicked() {
  QListWidgetItem *delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
  if (ui->listWidget->count() < 2) {
    ui->combinedOutputEdit->setEnabled(false);
  }
}

void topmgWindow::lockDialog() {
  ui->databaseFileButton->setEnabled(false);
  ui->modFileButton->setEnabled(false);
  ui->databaseFileEdit->setEnabled(false);
  ui->combinedOutputEdit->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  ui->fixedModFileEdit->setEnabled(false);
  ui->fixedModFileButton->setEnabled(false);
  ui->errorToleranceEdit->setEnabled(false);
  ui->maxModEdit->setEnabled(false);
  ui->cutoffSpectralValueEdit->setEnabled(false);
  ui->cutoffProteoformValueEdit->setEnabled(false);
  ui->modFileEdit->setEnabled(false);
  ui->threadNumberEdit->setEnabled(false);
  ui->skipListEdit->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  ui->activationComboBox->setEnabled(false);
  ui->cutoffSpectralTypeComboBox->setEnabled(false);
  ui->cutoffProteoformTypeComboBox->setEnabled(false);
  ui->numModComboBox->setEnabled(false);
  ui->numUnknownShiftComboBox->setEnabled(false);
  ui->NONECheckBox->setEnabled(false);
  ui->NMECheckBox->setEnabled(false);
  ui->NMEACCheckBox->setEnabled(false);
  ui->MACCheckBox->setEnabled(false);
  ui->decoyCheckBox->setEnabled(false);
  ui->topfdFeatureCheckBox->setEnabled(false);
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->outputButton->setEnabled(false);
  ui->addButton->setEnabled(false);
  ui->delButton->setEnabled(false);
}

void topmgWindow::unlockDialog() {
  ui->databaseFileButton->setEnabled(true);
  ui->modFileButton->setEnabled(true);
  ui->databaseFileEdit->setEnabled(true);
  if (ui->listWidget->count() > 1) {
    ui->combinedOutputEdit->setEnabled(true);
  }
  ui->fixedModFileEdit->setEnabled(true);
  ui->errorToleranceEdit->setEnabled(true);
  ui->maxModEdit->setEnabled(true);
  ui->cutoffSpectralValueEdit->setEnabled(true);
  ui->cutoffProteoformValueEdit->setEnabled(true);
  ui->modFileEdit->setEnabled(true);
  ui->threadNumberEdit->setEnabled(true);
  ui->skipListEdit->setEnabled(true);
  ui->fixedModComboBox->setEnabled(true);
  on_fixedModComboBox_currentIndexChanged(ui->fixedModComboBox->currentIndex());
  ui->activationComboBox->setEnabled(true);
  ui->cutoffSpectralTypeComboBox->setEnabled(true);
  ui->cutoffProteoformTypeComboBox->setEnabled(true);
  ui->numModComboBox->setEnabled(true);
  ui->numUnknownShiftComboBox->setEnabled(true);
  ui->NONECheckBox->setEnabled(true);
  ui->NMECheckBox->setEnabled(true);
  ui->NMEACCheckBox->setEnabled(true);
  ui->MACCheckBox->setEnabled(true);
  ui->decoyCheckBox->setEnabled(true);
  ui->topfdFeatureCheckBox->setEnabled(true);
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);
  ui->addButton->setEnabled(true);
  ui->delButton->setEnabled(true);
}

bool topmgWindow::checkError() {
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
                         tr("The output folder path is too long!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->fixedModFileEdit->text().isEmpty() && ui->fixedModComboBox->currentIndex() == 3) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Please select a fixed modification file!"),
                         QMessageBox::Yes);
    return true;
  }

  if (ui->modFileEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Please select a variable modification file!"),
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
  if (ui->threadNumberEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Thread number is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  return false;
}


void topmgWindow::updateMsg(std::string msg) {
  showInfo = msg.c_str();
  ui->outputTextBrowser->setText(showInfo);
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
  QString info = msg.c_str();
}

void topmgWindow::showArguments() {
  QMessageBox::warning(0, "Arguments", ("executiveDir:" + arguments_["executiveDir"] +
                                        "\nresourceDir:" + arguments_["resourceDir"] +
                                        "\noriDatabaseFileName:" + arguments_["oriDatabaseFileName"] +
                                        "\ndatabaseFileName:" + arguments_["databaseFileName"] +
                                        "\ndatabaseBlockSize:" + arguments_["databaseBlockSize"] +
                                        "\nspectrumFileName:" + arguments_["spectrumFileName"] +
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
                                        "\nkeepTempFiles:" + arguments_["keepTempFiles"] +
                                        "\ngroupSpectrumNumber:" + arguments_["groupSpectrumNumber"] +
                                        "\nfilteringResultNumber:" + arguments_["filteringResultNumber"] +
                                        "\nvarModFileName:" + arguments_["varModFileName"] +
                                        "\nthreadNumber:" + arguments_["threadNumber"] +
                                        "\nuseFeatureFile:" + arguments_["useFeatureFile"] +
                                        "\nskipList:" + arguments_["skipList"] +
                                        "\nproteoGraphGap:" + arguments_["proteoGraphGap"] +
                                        "\nuseAsfDiag:" + arguments_["useAsfDiag"] +
                                        "\nvarPtmNumber:" + arguments_["varPtmNumber"] +
                                        "\nvarPtmNumInGap:" + arguments_["varPtmNumInGap"]).c_str(), QMessageBox::Yes);
}

void topmgWindow::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait)
    QCoreApplication::processEvents();
}

void topmgWindow::on_exitButton_clicked() {
  close();
}

void topmgWindow::on_fixedModComboBox_currentIndexChanged(int index) {
  if (index == 3) {
    ui->fixedModFileEdit->setEnabled(true);
    ui->fixedModFileButton->setEnabled(true);
  } else {
    ui->fixedModFileEdit->setEnabled(false);
    ui->fixedModFileButton->setEnabled(false);
  }
}

bool topmgWindow::nterminalerror() {
  if (ui->NONECheckBox->isChecked() || ui->NMECheckBox->isChecked() || ui->NMEACCheckBox->isChecked() || ui->MACCheckBox->isChecked()) {
    return false;
  } else {
    QMessageBox::warning(this, tr("Warning"),
                         tr("At least one N-terminal form should be selected!"),
                         QMessageBox::Yes);
    return true;
  }
}

void topmgWindow::on_NONECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NONECheckBox->setChecked(true);
  }
}

void topmgWindow::on_NMECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMECheckBox->setChecked(true);
  }
}

void topmgWindow::on_NMEACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMEACCheckBox->setChecked(true);
  }
}

void topmgWindow::on_MACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->MACCheckBox->setChecked(true);
  }
}

void topmgWindow::on_cutoffSpectralTypeComboBox_currentIndexChanged(int index) {
  if (index == 1 && !ui->decoyCheckBox->isChecked()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an FDR cutoff, the checkbox \"decoy database\" should be checked!"),
                         QMessageBox::Yes);
    ui->cutoffSpectralTypeComboBox->setCurrentIndex(0);
  }
}

void topmgWindow::on_cutoffProteoformTypeComboBox_currentIndexChanged(int index) {
  if (index == 1 && !ui->decoyCheckBox->isChecked()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an FDR cutoff, the checkbox \"decoy database\" should be checked!"),
                         QMessageBox::Yes);
    ui->cutoffProteoformTypeComboBox->setCurrentIndex(0);
  }
}

void topmgWindow::on_decoyCheckBox_clicked(bool checked) {
  if (!checked && (ui->cutoffSpectralTypeComboBox->currentIndex() > 0 || ui->cutoffProteoformTypeComboBox->currentIndex() > 0)) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Because an FDR cutoff is selected, the checkbox \"decoy database\" cannot be unchecked."),
                         QMessageBox::Yes);
    ui->decoyCheckBox->setChecked(true);
  }
}

void topmgWindow::closeEvent(QCloseEvent *event) {
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

bool topmgWindow::continueToClose() {
  if (QMessageBox::question(this,
                            tr("Quit"),
                            tr("TopMG is still running. Are you sure you want to quit?"),
                            QMessageBox::Yes | QMessageBox::No,
                            QMessageBox::No)
      == QMessageBox::Yes) {
    return true;
  } else {
    return false;
  }
}
