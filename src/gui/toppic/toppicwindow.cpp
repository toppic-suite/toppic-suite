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

#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QCloseEvent>
#include <QCoreApplication>
#include <QToolTip>
#include <QDesktopServices>
#include <QScrollBar>
#include <QDebug>

#include "common/util/version.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/file_util.hpp"
#include "console/toppic_argument.hpp"

#include "gui/util/command.hpp"
#include "gui/util/gui_message.hpp"

#include "gui/toppic/ui_toppicwindow.h"
#include "gui/toppic/toppicwindow.hpp"

ToppicWindow::ToppicWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ToppicWindow) {
      ui->setupUi(this);
      std::string title = "TopPIC v." + toppic::Version::getVersion();
      QString qstr = QString::fromStdString(title);
      this->setWindowTitle(qstr);
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
      ui->formErrorToleranceEdit->setValidator(new QDoubleValidator(0, 2147483647, 4, this));
      QRegExp rx3("^-?\\d{1,8}\\.\\d{0,2}$");
      QRegExpValidator *validator3 = new QRegExpValidator(rx3, this);
      ui->minModEdit->setValidator(validator3);

      QFont font;
      QFont outputFont;
      QFont tableFont;
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
      font.setFamily(QStringLiteral("Calibri"));
      tableFont.setFamily(QStringLiteral("Calibri"));
      outputFont.setFamily(QStringLiteral("Consolas"));
#else
      font.setFamily(QStringLiteral("Monospace"));
      tableFont.setFamily(QStringLiteral("Monospace"));
      outputFont.setFamily(QStringLiteral("Monospace"));
#endif
      font.setPixelSize(12);
      QApplication::setFont(font);
      outputFont.setPixelSize(12);
      ui->outputTextBrowser->setFont(outputFont);
      tableFont.setPointSize(9);
      ui->listWidget->setFont(tableFont);

      setToolTip("");
      setToolTipDuration(100);

      on_clearButton_clicked();
      on_defaultButton_clicked();
      ui->tabWidget->setCurrentIndex(0);
    }

ToppicWindow::~ToppicWindow() {
  if(process_.state()!=QProcess::NotRunning) {
    process_.kill();
  }
  delete ui;
}

void ToppicWindow::on_clearButton_clicked() {
  ui->databaseFileEdit->clear();
  ui->listWidget->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
  ui->combinedOutputEdit->setText("");
  ui->combinedOutputEdit->setEnabled(false);
  ui->outputButton->setEnabled(false);
}

void ToppicWindow::on_defaultButton_clicked() {
  arguments_ = toppic::ToppicArgument::initArguments();
  
  ui->combinedOutputEdit->setText("");
  ui->fixedModFileEdit->clear();
  ui->errorToleranceEdit->setText(QString::fromStdString(arguments_["massErrorTolerance"])); 
  ui->formErrorToleranceEdit->setText(QString::fromStdString(arguments_["proteoformErrorTolerance"]));
  ui->varPtmNumEdit->setText(QString::fromStdString(arguments_["variablePtmNum"])); 
  ui->varPtmFileEdit->clear();
  ui->approxSpectraCheckBox->setChecked(false);
  ui->minModEdit->setText(QString::fromStdString(arguments_["minShiftMass"]));
  ui->maxModEdit->setText(QString::fromStdString(arguments_["maxShiftMass"]));
  ui->cutoffSpectralValueEdit->setText(QString::fromStdString(arguments_["cutoffSpectralValue"]));
  ui->cutoffProteoformValueEdit->setText(QString::fromStdString(arguments_["cutoffProteoformValue"]));
  ui->numCombinedEdit->setText(QString::fromStdString(arguments_["groupSpectrumNumber"]));
  ui->miscoreThresholdEdit->setText(QString::fromStdString(arguments_["localThreshold"]));
  ui->threadNumberEdit->setText(QString::fromStdString(arguments_["threadNumber"]));

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
  ui->geneHTMLCheckBox->setChecked(true);
  ui->keepDecoyCheckBox->setChecked(false);
  ui->keepTempCheckBox->setChecked(false);
}

void ToppicWindow::updatedir(QString s) {
  if (!s.isEmpty()) {
    //lastDir_ = s;
    lastDir_ = "";
  }
}

// single file
void ToppicWindow::ToppicWindow::on_databaseFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a protein database file",
      lastDir_,
      "Database files (*.fasta *.fa)");
  updatedir(s);
  ui->databaseFileEdit->setText(s);
}

void ToppicWindow::on_fixedModFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a fixed modification file",
      lastDir_,
      "Modification files (*.txt);;All files (*.*)");
  updatedir(s);
  ui->fixedModFileEdit->setText(s);
}

void ToppicWindow::on_modFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a modification file for PTM localization",
      lastDir_,
      "Modification files (*.txt);;All files (*.*)");
  updatedir(s);
  ui->modFileEdit->setText(s);
}

void ToppicWindow::on_varPtmFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a modification file for variable PTMs",
      lastDir_,
      "Modification files (*.txt);;All files (*.*)");
  updatedir(s);
  ui->varPtmFileEdit->setText(s);
}

void ToppicWindow::on_startButton_clicked() {
  lockDialog();

  std::map<std::string, std::string> argument = this->getArguments();
  std::vector<std::string> spec_file_lst = this->getSpecFileList();

  std::string cmd = toppic::command::geneToppicCommand(argument, spec_file_lst);
  QString q_cmd = QString::fromStdString(cmd);
  q_cmd = q_cmd.trimmed();
  QStringList cmd_list = q_cmd.split(" ");
  QString prog = cmd_list[0];
  cmd_list.removeFirst();
  //qDebug() << q_cmd;

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

void ToppicWindow::on_outputButton_clicked() {
  std::vector<std::string> spec_file_lst = this->getSpecFileList();
  if (spec_file_lst.size() > 0) {
    std::string dir = toppic::file_util::directory(spec_file_lst[0]);
    QString outPath = dir.c_str();
    QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
  }
}

std::map<std::string, std::string> ToppicWindow::getArguments() {
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
  arguments_["shiftNumber"] = ui->numModComboBox->currentText().toStdString();
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

  arguments_["minShiftMass"] = ui->minModEdit->text().toStdString();
  arguments_["maxShiftMass"] = ui->maxModEdit->text().toStdString();

  if (ui->approxSpectraCheckBox->isChecked()) {
    arguments_["useApproxSpectra"] = "true";
  } else {
    arguments_["useApproxSpectra"] = "false";
  }
  arguments_["variablePtmNum"] = ui->varPtmNumEdit->text().toStdString();
  arguments_["variablePtmFileName"] = ui->varPtmFileEdit->text().trimmed().toStdString();

  if (ui->lookupTableCheckBox->isChecked()) {
    arguments_["useLookupTable"] = "true";
  } else {
    arguments_["useLookupTable"] = "false";
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
  arguments_["localThreshold"] = ui->miscoreThresholdEdit->text().toStdString();
  arguments_["groupSpectrumNumber"] = ui->numCombinedEdit->text().toStdString();
  arguments_["localPtmFileName"] = ui->modFileEdit->text().toStdString();
  arguments_["threadNumber"] = ui->threadNumberEdit->text().toStdString();
  if (ui->topfdFeatureCheckBox->isChecked()) {
    arguments_["useFeatureFile"] = "false";
  } else {
    arguments_["useFeatureFile"] = "true";
  }
  if (ui->geneHTMLCheckBox->isChecked()) {
    arguments_["geneHTMLFolder"] = "true";
  } else {
    arguments_["geneHTMLFolder"] = "false";
  }
  //showArguments();
  return arguments_;
}

std::vector<std::string> ToppicWindow::getSpecFileList() {
  spec_file_lst_.clear();
  for (int i = 0; i < ui->listWidget->count(); i++) {
    spec_file_lst_.push_back(ui->listWidget->item(i)->text().toStdString());
  }
  return spec_file_lst_;
}

void ToppicWindow::on_addButton_clicked() {
  QString filter = "Spectrum files (*ms2.msalign)";
  QStringList spfiles = QFileDialog::getOpenFileNames(
      this,
      "Select deconvoluted spectrum files",
      lastDir_,
      filter
      );

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

bool ToppicWindow::ableToAdd(QString spfile) {
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

void ToppicWindow::on_delButton_clicked() {
  QListWidgetItem *delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
  if (ui->listWidget->count() < 2) {
    ui->combinedOutputEdit->setEnabled(false);
  }
}

void ToppicWindow::lockDialog() {
  ui->databaseFileButton->setEnabled(false);
  ui->modFileButton->setEnabled(false);
  ui->databaseFileEdit->setEnabled(false);
  ui->combinedOutputEdit->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  ui->fixedModFileEdit->setEnabled(false);
  ui->fixedModFileButton->setEnabled(false);
  ui->errorToleranceEdit->setEnabled(false);
  ui->formErrorToleranceEdit->setEnabled(false);
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
  ui->geneHTMLCheckBox->setEnabled(false);
  ui->keepDecoyCheckBox->setEnabled(false);
  ui->keepTempCheckBox->setEnabled(false);
}

void ToppicWindow::unlockDialog() {
  ui->databaseFileButton->setEnabled(true);
  ui->modFileButton->setEnabled(true);
  ui->databaseFileEdit->setEnabled(true);
  if (ui->listWidget->count() > 1) {
    ui->combinedOutputEdit->setEnabled(true);
  }
  ui->fixedModFileEdit->setEnabled(true);
  ui->errorToleranceEdit->setEnabled(true);
  ui->formErrorToleranceEdit->setEnabled(true);
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
  ui->geneHTMLCheckBox->setEnabled(true);
  ui->keepDecoyCheckBox->setEnabled(true);
  ui->keepTempCheckBox->setEnabled(true);
}

bool ToppicWindow::checkError() {
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
                         tr("Mass error tolerance is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->formErrorToleranceEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Proteoform error tolerance is empty!"),
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

void ToppicWindow::updateMsg(std::string msg) {
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

void ToppicWindow::showArguments() {
  QMessageBox::warning(0, "Arguments", ("executiveDir:" + arguments_["executiveDir"] +
                                        "\nresourceDir:" + arguments_["resourceDir"] +
                                        "\noriDatabaseFileName:" + arguments_["oriDatabaseFileName"] +
                                        "\ndatabaseFileName:" + arguments_["databaseFileName"] +
                                        "\ndatabaseBlockSize:" + arguments_["databaseBlockSize"] +
                                        "\nmaxFragmentLength:" + arguments_["maxFragmentLength"] +
                                        "\ncombinedOutputName:" + arguments_["combinedOutputName"] +
                                        "\nactivation:" + arguments_["activation"] +
                                        "\nsearchType:" + arguments_["searchType"] +
                                        "\nfixedMod:" + arguments_["fixedMod"] +
                                        "\nptmNumber:" + arguments_["ptmNumber"] +
                                        "\nmassErrorTolerance:" + arguments_["massErrorTolerance"] +
                                        "\ncutoffSpectralType:" + arguments_["cutoffSpectralType"] +
                                        "\ncutoffSpectralValue:" + arguments_["cutoffSpectralValue"] +
                                        "\ncutoffProteoformType:" + arguments_["cutoffProteoformType"] +
                                        "\ncutoffProteoformValue:" + arguments_["cutoffProteoformValue"] +
                                        "\nallowProtMod:" + arguments_["allowProtMod"] +
                                        "\nnumOfTopPrsms:" + arguments_["numOfTopPrsms"] +
                                        "\nminShiftMass:" + arguments_["minShiftMass"] +
                                        "\nmaxShiftMass:" + arguments_["maxShiftMass"] +
                                        "\nuseLookupTable:" + arguments_["useLookupTable"] +
                                        "\nkeepTempFiles:" + arguments_["keepTempFiles"] +
                                        "\nlocalThreshold:" + arguments_["localThreshold"] +
                                        "\ngroupSpectrumNumber:" + arguments_["groupSpectrumNumber"] +
                                        "\nfilteringResultNumber:" + arguments_["filteringResultNumber"] +
                                        "\nlocalPtmFileName:" + arguments_["localPtmFileName"] +
                                        "\nthreadNumber:" + arguments_["threadNumber"] +
                                        "\nuseFeatureFile:" + arguments_["useFeatureFile"] +
                                        "\nskipList:" + arguments_["skipList"]).c_str(), QMessageBox::Yes);
}

void ToppicWindow::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait)
    QCoreApplication::processEvents();
}

void ToppicWindow::on_exitButton_clicked() {
  close();
}

void ToppicWindow::on_fixedModComboBox_currentIndexChanged(int index) {
  if (index == 3) {
    ui->fixedModFileEdit->setEnabled(true);
    ui->fixedModFileButton->setEnabled(true);
  } else {
    ui->fixedModFileEdit->setEnabled(false);
    ui->fixedModFileButton->setEnabled(false);
  }
}

void ToppicWindow::on_numModComboBox_currentIndexChanged(int index) {
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

void ToppicWindow::on_errorToleranceEdit_textChanged(QString string) {
  QString currentText = ui->errorToleranceEdit->text();
  if (ui->lookupTableCheckBox->isChecked() && currentText != "5" && currentText != "10" && currentText != "15" && currentText != "1") {
    QMessageBox::warning(this, tr("Warning"),
                         tr("When the checkbox \"Lookup table for E-value computation\" is checked, only three error tolerance values 5, 10, and 15 ppm can be used!"),
                         QMessageBox::Yes);
    ui->errorToleranceEdit->setText("15");
  }
}

void ToppicWindow::on_lookupTableCheckBox_clicked(bool checked) {
  QString currentText = ui->errorToleranceEdit->text();
  if (checked && currentText != "5" && currentText != "10" && currentText != "15") {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an error tolerance other than 5, 10, and 15 ppm, the checkbox \"Lookup table for E-value computation\" should not be checked!"),
                         QMessageBox::Yes);
    ui->lookupTableCheckBox->setChecked(true);
  }
}

bool ToppicWindow::nterminalerror() {
  if (ui->NONECheckBox->isChecked() || ui->NMECheckBox->isChecked() || ui->NMEACCheckBox->isChecked() || ui->MACCheckBox->isChecked()) {
    return false;
  } else {
    QMessageBox::warning(this, tr("Warning"),
                         tr("At least one N-terminal form should be selected!"),
                         QMessageBox::Yes);
    return true;
  }
}

void ToppicWindow::on_NONECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NONECheckBox->setChecked(true);
  }
}

void ToppicWindow::on_NMECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMECheckBox->setChecked(true);
  }
}

void ToppicWindow::on_NMEACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMEACCheckBox->setChecked(true);
  }
}

void ToppicWindow::on_MACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->MACCheckBox->setChecked(true);
  }
}

void ToppicWindow::on_cutoffSpectralTypeComboBox_currentIndexChanged(int index) {
  if (index == 1 && !ui->decoyCheckBox->isChecked()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an FDR cutoff, the checkbox \"decoy database\" should be checked!"),
                         QMessageBox::Yes);
    ui->cutoffSpectralTypeComboBox->setCurrentIndex(0);
  }
}

void ToppicWindow::on_cutoffProteoformTypeComboBox_currentIndexChanged(int index) {
  if (index == 1 && !ui->decoyCheckBox->isChecked()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an FDR cutoff, the checkbox \"decoy database\" should be checked!"),
                         QMessageBox::Yes);
    ui->cutoffProteoformTypeComboBox->setCurrentIndex(0);
  }
}

void ToppicWindow::on_decoyCheckBox_clicked(bool checked) {
  if (!checked && (ui->cutoffSpectralTypeComboBox->currentIndex() > 0 || ui->cutoffProteoformTypeComboBox->currentIndex() > 0)) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Because an FDR cutoff is selected, the checkbox \"decoy database\" cannot be unchecked."),
                         QMessageBox::Yes);
    ui->decoyCheckBox->setChecked(true);
  }
}

void ToppicWindow::closeEvent(QCloseEvent *event) {
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
bool ToppicWindow::continueToClose() {
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

bool ToppicWindow::event(QEvent *event) {
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
