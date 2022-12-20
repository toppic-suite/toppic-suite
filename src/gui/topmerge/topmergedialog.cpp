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

#include "common/util/file_util.hpp"
#include "common/base/base_data.hpp"
#include "common/util/version.hpp"

#include "gui/topmerge/topmergedialog.h"
#include "gui/topmerge/ui_topmergedialog.h"
#include "gui/topmerge/threadtopmerge.h"


TopMergeDialog::TopMergeDialog(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::TopMergeDialog) {
      initArguments();
      ui->setupUi(this);
      std::string title = "TopMerge v." + toppic::Version::getVersion();
      QString qstr = QString::fromStdString(title);
      this->setWindowTitle(qstr);
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
      thread_ = new ThreadTopMerge(this);
      showInfo = "";
      TopMergeDialog::on_defaultButton_clicked();
    }

TopMergeDialog::~TopMergeDialog() {
  thread_->terminate();
  thread_->wait();
  delete ui;
}

void TopMergeDialog::on_databaseFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a protein database file",
      lastDir_,
      "Database files(*.fasta *.fa)");
  updatedir(s);
  ui->databaseFileEdit->setText(s);
}

void TopMergeDialog::closeEvent(QCloseEvent *event) {
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

void TopMergeDialog::initArguments() {
  arguments_["oriDatabaseFileName"]="";
  arguments_["databaseFileName"] = "";
  arguments_["combinedOutputName"] = "";
  arguments_["proteoformErrorTolerance"] = "1.2";
  arguments_["cutoffSpectralType"] = "EVALUE";
  arguments_["cutoffSpectralValue"] = "0.01";
  arguments_["cutoffProteoformType"] = "EVALUE";
  arguments_["cutoffProteoformValue"] = "0.01";
  arguments_["maxPtmMass"] = "500";
  arguments_["minPtmMass"] = "-500";
  arguments_["keepTempFiles"] = "false";
  arguments_["keepDecoyResults"] = "false";
  arguments_["localThreshold"] = "0.15";
  arguments_["residueModFileName"] = "";
  arguments_["useFeatureFile"] = "true";
  arguments_["geneHTMLFolder"] = "";
  arguments_["searchType"] = "TARGET";
  arguments_["massErrorTolerance"] = "15";
  arguments_["executiveDir"] = ".";
  arguments_["resourceDir"] = "";
  arguments_["activation"] = "FILE";
  arguments_["allowProtMod"] = "NONE,NME,NME_ACETYLATION,M_ACETYLATION";
  arguments_["fixedMod"] = "";
}

void TopMergeDialog::on_clearButton_clicked() {
  ui->databaseFileEdit->clear();
  ui->listWidget->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the data.");
  ui->combinedOutputEdit->setText("combined");
  ui->outputButton->setEnabled(true);
}

void TopMergeDialog::on_defaultButton_clicked() {
  ui->combinedOutputEdit->setText("combined");
  ui->fixedModFileEdit->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the data.");
  ui->errorToleranceEdit->setText("15");
  ui->formErrorToleranceEdit->setText("1.2");
  ui->maxModEdit->setText("500");
  ui->minModEdit->setText("-500");
  ui->cutoffSpectralValueEdit->setText("0.01");
  ui->cutoffProteoformValueEdit->setText("0.01");
  ui->miscoreThresholdEdit->setText("0.15");
  ui->fixedModComboBox->setCurrentIndex(0);
  on_fixedModComboBox_currentIndexChanged(0);
  ui->activationComboBox->setCurrentIndex(0);
  ui->cutoffSpectralTypeComboBox->setCurrentIndex(0);
  ui->cutoffProteoformTypeComboBox->setCurrentIndex(0);
  ui->topfdFeatureCheckBox->setChecked(false);
  ui->geneHTMLCheckBox->setChecked(true);
  ui->decoyCheckBox->setChecked(false);
  ui->keepDecoyCheckBox->setChecked(false);
  ui->keepTempCheckBox->setChecked(false);
  ui->NONECheckBox->setChecked(true);
  ui->NMECheckBox->setChecked(true);
  ui->NMEACCheckBox->setChecked(true);
  ui->MACCheckBox->setChecked(true);
}

void TopMergeDialog::updatedir(QString s) {
  if (!s.isEmpty()) {
    lastDir_ = s;
  }
}

void TopMergeDialog::on_startButton_clicked() {
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
  bool finish = false;

  while (true) {
    if (thread_->isFinished()) {
      finish = true;
    }
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
    if (finish) {
      break;
    }
    sleep(100);
  }
  unlockDialog();
  showInfo = "";
  thread_->exit();
  std::cout.rdbuf(oldbuf);
  
}

void TopMergeDialog::on_exitButton_clicked() {
  close();
}

bool TopMergeDialog::continueToClose() {
  if (QMessageBox::question(this,
                            tr("Quit"),
                            tr("TopMerge is still running. Are you sure you want to quit?"),
                            QMessageBox::Yes | QMessageBox::No,
                            QMessageBox::No)
      == QMessageBox::Yes) {
    return true;
  } else {
    return false;
  }
}

void TopMergeDialog::on_fixedModFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a fixed modification file",
      lastDir_,
      "Modification files(*.txt);;All files(*.*)");
  updatedir(s);
  ui->fixedModFileEdit->setText(s);
}

void TopMergeDialog::on_fixedModComboBox_currentIndexChanged(int index) {
  if (index == 3) {
    ui->fixedModFileEdit->setEnabled(true);
    ui->fixedModFileButton->setEnabled(true);
  } else {
    ui->fixedModFileEdit->setEnabled(false);
    ui->fixedModFileButton->setEnabled(false);
  }
}
bool TopMergeDialog::nterminalerror() {
  if (ui->NONECheckBox->isChecked() || ui->NMECheckBox->isChecked() || ui->NMEACCheckBox->isChecked() || ui->MACCheckBox->isChecked()) {
    return false;
  } else {
    QMessageBox::warning(this, tr("Warning"),
                         tr("At least one N-terminal form should be selected!"),
                         QMessageBox::Yes);
    return true;
  }
}

void TopMergeDialog::on_NONECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NONECheckBox->setChecked(true);
  }
}

void TopMergeDialog::on_NMECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMECheckBox->setChecked(true);
  }
}

void TopMergeDialog::on_NMEACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMEACCheckBox->setChecked(true);
  }
}

void TopMergeDialog::on_MACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->MACCheckBox->setChecked(true);
  }
}

void TopMergeDialog::on_modFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a modification file for PTM localization",
      lastDir_,
      "Modification files(*.txt);;All files(*.*)");
  updatedir(s);
  ui->modFileEdit->setText(s);
}

void TopMergeDialog::on_outputButton_clicked() {
  std::string db_file_name = ui->databaseFileEdit->text().toStdString();

  std::string dir = toppic::file_util::directory(db_file_name);
  QString outPath = dir.c_str();

  QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
}

std::map<std::string, std::string> TopMergeDialog::getArguments() {
  QString path = QCoreApplication::applicationFilePath();
  arguments_["executiveDir"] = toppic::file_util::getExecutiveDir(path.toStdString());
  if (toppic::file_util::checkSpace(arguments_["executiveDir"])) {
    ui->outputTextBrowser->setText("Current directory " + QString::fromStdString(arguments_["executiveDir"]) + " contains space and will cause errors in the program!");
  }
  arguments_["resourceDir"] = toppic::file_util::getResourceDir(arguments_["executiveDir"]);
  arguments_["oriDatabaseFileName"] = ui->databaseFileEdit->text().toStdString();
  arguments_["combinedOutputName"] = ui->combinedOutputEdit->text().trimmed().toStdString();
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
  else if (arguments_["fixedMod"] == "Carbamidomethylation on cysteine") {
    arguments_["fixedMod"] = "C57";
  }
  else if (arguments_["fixedMod"] == "Carboxymethylation on cysteine") {
    arguments_["fixedMod"] = "C58";
  }
  if (ui->fixedModComboBox->currentIndex() == 3) {
    arguments_["fixedMod"] = ui->fixedModFileEdit->text().toStdString();
  }
  arguments_["massErrorTolerance"] = ui->errorToleranceEdit->text().toStdString();
  arguments_["proteoformErrorTolerance"] = ui->formErrorToleranceEdit->text().toStdString();
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
  arguments_["maxPtmMass"] = ui->maxModEdit->text().toStdString();
  arguments_["minPtmMass"] = ui->minModEdit->text().toStdString();
  arguments_["keepTempFiles"] = "false";   // default
  arguments_["localThreshold"] = ui->miscoreThresholdEdit->text().toStdString();
  arguments_["residueModFileName"] = ui->modFileEdit->text().toStdString();
  if (ui->topfdFeatureCheckBox->isChecked()) {
    arguments_["useFeatureFile"] = "false";
  } else {
    arguments_["useFeatureFile"] = "true";
  }
  if (ui->keepTempCheckBox->isChecked()) {
    arguments_["keepTempFiles"] = "true";
  } else {
    arguments_["keepTempFiles"] = "false";
  }
  if (ui->keepDecoyCheckBox->isChecked()) {
    arguments_["keepDecoyResults"] = "true";
  } else {
    arguments_["keepDecoyResults"] = "false";
  }
  if (ui->geneHTMLCheckBox->isChecked()) {
    arguments_["geneHTMLFolder"] = "true";
  } else {
    arguments_["geneHTMLFolder"] = "false";
  }
  return arguments_;
}
std::vector<std::string> TopMergeDialog::getSpecFileList() {
  spec_file_lst_.clear();
  for (int i = 0; i < ui->listWidget->count(); i++) {
    spec_file_lst_.push_back(ui->listWidget->item(i)->text().toStdString());
  }
  return spec_file_lst_;
}
void TopMergeDialog::on_addButton_clicked() {
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

bool TopMergeDialog::ableToAdd(QString idfile) {
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
void TopMergeDialog::on_delButton_clicked() {
  QListWidgetItem *delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
}
void TopMergeDialog::lockDialog() {
  ui->combinedOutputEdit->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  ui->fixedModFileEdit->setEnabled(false);
  ui->fixedModFileButton->setEnabled(false);
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->outputButton->setEnabled(false);
  ui->errorToleranceEdit->setEnabled(false);
  ui->formErrorToleranceEdit->setEnabled(false);  
  ui->maxModEdit->setEnabled(false);  
  ui->minModEdit->setEnabled(false);  
  ui->cutoffSpectralValueEdit->setEnabled(false);  
  ui->cutoffProteoformValueEdit->setEnabled(false);  
  ui->modFileEdit->setEnabled(false);
  ui->miscoreThresholdEdit->setEnabled(false);  
  ui->cutoffSpectralTypeComboBox->setEnabled(false);  
  ui->cutoffProteoformTypeComboBox->setEnabled(false);  
  ui->databaseFileButton->setEnabled(false);
  ui->databaseFileEdit->setEnabled(false);
  ui->decoyCheckBox->setEnabled(false);
  ui->topfdFeatureCheckBox->setEnabled(false);
  ui->addButton->setEnabled(false);
  ui->delButton->setEnabled(false);
  ui->geneHTMLCheckBox->setEnabled(false);
  ui->keepDecoyCheckBox->setEnabled(false);
  ui->keepTempCheckBox->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  ui->activationComboBox->setEnabled(false);
  ui->NONECheckBox->setEnabled(false);
  ui->NMECheckBox->setEnabled(false);
  ui->NMEACCheckBox->setEnabled(false);
  ui->MACCheckBox->setEnabled(false);
}

void TopMergeDialog::unlockDialog() {
  ui->combinedOutputEdit->setEnabled(true);
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);

  ui->databaseFileButton->setEnabled(true);
  ui->databaseFileEdit->setEnabled(true);
  ui->fixedModFileEdit->setEnabled(true);
  ui->errorToleranceEdit->setEnabled(true);  
  ui->formErrorToleranceEdit->setEnabled(true);  
  ui->maxModEdit->setEnabled(true);  
  ui->minModEdit->setEnabled(true);  
  ui->fixedModComboBox->setEnabled(true);
  on_fixedModComboBox_currentIndexChanged(ui->fixedModComboBox->currentIndex());
  ui->activationComboBox->setEnabled(true);
  ui->cutoffSpectralValueEdit->setEnabled(true);  
  ui->cutoffProteoformValueEdit->setEnabled(true);  
  ui->miscoreThresholdEdit->setEnabled(true);  
  ui->cutoffSpectralTypeComboBox->setEnabled(true);  
  ui->cutoffProteoformTypeComboBox->setEnabled(true);  
  ui->NONECheckBox->setEnabled(true);
  ui->NMECheckBox->setEnabled(true);
  ui->NMEACCheckBox->setEnabled(true);
  ui->MACCheckBox->setEnabled(true);
  ui->decoyCheckBox->setEnabled(true);
  ui->topfdFeatureCheckBox->setEnabled(true);
  ui->geneHTMLCheckBox->setEnabled(true);
  ui->keepDecoyCheckBox->setEnabled(true);
  ui->keepTempCheckBox->setEnabled(true);
  ui->addButton->setEnabled(true);
  ui->delButton->setEnabled(true);
}

bool TopMergeDialog::checkError() {
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
  if (ui->listWidget->count() == 0) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Spectrum files are not selected!"),
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
                         tr("Mass error tolerance is empty!"),
                         QMessageBox::Yes);
    return true;
  }

  return false;
}

void TopMergeDialog::updateMsg(std::string msg) {
  showInfo = msg.c_str();
  ui->outputTextBrowser->setText(showInfo);
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
}

void TopMergeDialog::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait) {
    QCoreApplication::processEvents();
  }
}