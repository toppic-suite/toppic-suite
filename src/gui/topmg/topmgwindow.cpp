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

#include "common/util/version.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/file_util.hpp"

#include "console/topmg_argument.hpp"

#include "gui/util/command.hpp"
#include "gui/util/gui_message.hpp"

#include "gui/topmg/ui_topmgwindow.h"
#include "gui/topmg/topmgwindow.hpp"

TopmgWindow::TopmgWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::TopmgWindow) {
      ui->setupUi(this);
      std::string title = "TopMG v." + toppic::Version::getVersion();
      QString qstr = QString::fromStdString(title);
      this->setWindowTitle(qstr);
      lastDir_ = ".";
      QRegExp rx1("^\\d{1,8}\\.\\d{0,2}$");
      QRegExpValidator *validator1 = new QRegExpValidator(rx1, this);
      ui->maxModEdit->setValidator(validator1);
      ui->cutoffSpectralValueEdit->setValidator(validator1);
      ui->cutoffProteoformValueEdit->setValidator(validator1);
      ui->threadNumberEdit->setValidator(new QIntValidator(0, 2147483647, this));
      ui->errorToleranceEdit->setValidator(new QIntValidator(0, 2147483647, this));
      ui->formErrorToleranceEdit->setValidator(new QDoubleValidator(0, 2147483647, 4, this));

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

TopmgWindow::~TopmgWindow() {
  if(process_.state()!=QProcess::NotRunning) {
    process_.kill();
  }
  delete ui;
}

void TopmgWindow::on_clearButton_clicked() {
  ui->databaseFileEdit->clear();
  ui->listWidget->clear();
  ui->combinedOutputEdit->clear();
  ui->combinedOutputEdit->setEnabled(false);
  ui->modFileEdit->setText("");
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum files.");
}

void TopmgWindow::on_defaultButton_clicked() {
  arguments_ = toppic::TopmgArgument::initArguments();
  
  ui->combinedOutputEdit->setText("");
  ui->modFileEdit->setText("");

  ui->errorToleranceEdit->setText(QString::fromStdString(arguments_["massErrorTolerance"])); 
  ui->formErrorToleranceEdit->setText(QString::fromStdString(arguments_["proteoformErrorTolerance"]));
  ui->maxModEdit->setText(QString::fromStdString(arguments_["maxShiftMass"]));
  ui->cutoffSpectralValueEdit->setText(QString::fromStdString(arguments_["cutoffSpectralValue"]));
  ui->cutoffProteoformValueEdit->setText(QString::fromStdString(arguments_["cutoffProteoformValue"]));
  ui->threadNumberEdit->setText(QString::fromStdString(arguments_["threadNumber"]));

  ui->fixedModComboBox->setCurrentIndex(0);
  ui->fixedModFileEdit->clear();
  ui->fixedModFileEdit->setEnabled(false);
  ui->fixedModFileButton->setEnabled(false);
  ui->activationComboBox->setCurrentIndex(0);
  ui->cutoffSpectralTypeComboBox->setCurrentIndex(0);
  ui->cutoffProteoformTypeComboBox->setCurrentIndex(0);
  // number of variable PTMs: 5
  ui->numModComboBox->setCurrentIndex(4);
  // number of unknown mass shifts: 0
  ui->numUnknownShiftComboBox->setCurrentIndex(0);
  ui->NONECheckBox->setChecked(true);
  ui->NMECheckBox->setChecked(true);
  ui->NMEACCheckBox->setChecked(true);
  ui->MACCheckBox->setChecked(true);
  ui->decoyCheckBox->setChecked(false);
  ui->topfdFeatureCheckBox->setChecked(false);
  ui->asfDiagCheckBox->setChecked(false);
  ui->geneHTMLCheckBox->setChecked(true);
  ui->wholeProteinCheckBox->setChecked(false);
  ui->maxGapLength->setText(QString::fromStdString(arguments_["proteoGraphGap"]));
  ui->maxVarPTMGap->setText(QString::fromStdString(arguments_["varPtmNumInGap"]));
  ui->keepDecoyCheckBox->setChecked(false);
  ui->keepTempCheckBox->setChecked(false);
}

void TopmgWindow::updatedir(QString s) {
  if (!s.isEmpty()) {
    //lastDir_ = s;
    lastDir_ = "";
  }
}

void TopmgWindow::TopmgWindow::on_databaseFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a protein database file",
      lastDir_,
      "Database files (*.fasta *.fa)");
  updatedir(s);
  ui->databaseFileEdit->setText(s);
}

void TopmgWindow::on_fixedModFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a fixed modification file",
      lastDir_,
      "Modification files (*.txt);;All files (*.*)");
  updatedir(s);
  ui->fixedModFileEdit->setText(s);
}

void TopmgWindow::on_modFileButton_clicked() {
  QString s = QFileDialog::getOpenFileName(
      this,
      "Select a modification file for variable PTMs",
      lastDir_,
      "Modification files (*.txt);;All files (*.*)");
  updatedir(s);
  ui->modFileEdit->setText(s);
}

void TopmgWindow::on_startButton_clicked() {
  lockDialog();
  std::map<std::string, std::string> argument = this->getArguments();
  std::vector<std::string> spec_file_lst = this->getSpecFileList();

  std::string cmd = toppic::command::geneTopmgCommand(argument, spec_file_lst);
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

void TopmgWindow::on_outputButton_clicked() {
  std::vector<std::string> spec_file_lst = this->getSpecFileList();
  if (spec_file_lst.size() > 0) {
    std::string dir = toppic::file_util::directory(spec_file_lst[0]);
    QString outPath = dir.c_str();
    QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
  }
}

std::map<std::string, std::string> TopmgWindow::getArguments() {
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
  arguments_["varPtmNumber"] = ui->numModComboBox->currentText().toStdString();
  arguments_["shiftNumber"] = ui->numUnknownShiftComboBox->currentText().toStdString();
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
  arguments_["maxShiftMass"] = ui->maxModEdit->text().toStdString();
  arguments_["proteoGraphGap"] = ui->maxGapLength->text().toStdString();
  arguments_["varPtmNumber"] = ui->numModComboBox->currentText().toStdString();
  arguments_["shiftNumber"] = ui->numUnknownShiftComboBox->currentText().toStdString();
  arguments_["varPtmNumInGap"] = ui->maxVarPTMGap->text().toStdString();
  arguments_["varModFileName"] = ui->modFileEdit->text().toStdString();
  arguments_["threadNumber"] = ui->threadNumberEdit->text().toStdString();
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
  if (ui->asfDiagCheckBox->isChecked()) {
    arguments_["useAsfDiag"] = "true";
  } else {
    arguments_["useAsfDiag"] = "false";
  }
  if (ui->geneHTMLCheckBox->isChecked()) {
    arguments_["geneHTMLFolder"] = "true";
  } else {
    arguments_["geneHTMLFolder"] = "false";
  }
  if (ui->wholeProteinCheckBox->isChecked()) {
    arguments_["wholeProteinOnly"] = "true";
  } else {
    arguments_["wholeProteinOnly"] = "false";
  }
  //showArguments();
  return arguments_;
}

std::vector<std::string> TopmgWindow::getSpecFileList() {
  spec_file_lst_.clear();
  for (int i = 0; i < ui->listWidget->count(); i++) {
    spec_file_lst_.push_back(ui->listWidget->item(i)->text().toStdString());
  }
  return spec_file_lst_;
}

void TopmgWindow::on_addButton_clicked() {
  QStringList spfiles = QFileDialog::getOpenFileNames(
      this,
      "Select deconvoluted spectrum files",
      lastDir_,
      "Spectrum files (*ms2.msalign)");
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

bool TopmgWindow::ableToAdd(QString spfile) {
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

void TopmgWindow::on_delButton_clicked() {
  QListWidgetItem *delItem = ui->listWidget->currentItem();
  ui->listWidget->removeItemWidget(delItem);
  delete delItem;
  if (ui->listWidget->count() < 2) {
    ui->combinedOutputEdit->setEnabled(false);
  }
}

void TopmgWindow::lockDialog() {
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
  ui->cutoffSpectralValueEdit->setEnabled(false);
  ui->cutoffProteoformValueEdit->setEnabled(false);
  ui->modFileEdit->setEnabled(false);
  ui->threadNumberEdit->setEnabled(false);
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
  ui->geneHTMLCheckBox->setEnabled(false);
  ui->wholeProteinCheckBox->setEnabled(false);
  ui->maxGapLength->setEnabled(false);
  ui->maxVarPTMGap->setEnabled(false);
  ui->keepDecoyCheckBox->setEnabled(false);
  ui->keepTempCheckBox->setEnabled(false);
  ui->asfDiagCheckBox->setEnabled(false);
}

void TopmgWindow::unlockDialog() {
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
  ui->cutoffSpectralValueEdit->setEnabled(true);
  ui->cutoffProteoformValueEdit->setEnabled(true);
  ui->modFileEdit->setEnabled(true);
  ui->threadNumberEdit->setEnabled(true);
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
  ui->geneHTMLCheckBox->setEnabled(true);
  ui->wholeProteinCheckBox->setEnabled(true);
  ui->maxGapLength->setEnabled(true);
  ui->maxVarPTMGap->setEnabled(true);
  ui->keepDecoyCheckBox->setEnabled(true);
  ui->keepTempCheckBox->setEnabled(true);
  ui->asfDiagCheckBox->setEnabled(true);
}

bool TopmgWindow::checkError() {
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
                         tr("Mass error tolerance is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->formErrorToleranceEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Prsm cluster error tolerance is empty!"),
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
  if (ui->threadNumberEdit->text().toInt() > toppic::mem_check::getMaxThreads("topmg")) {
    int max_thread = toppic::mem_check::getMaxThreads("topmg");
    QMessageBox::StandardButton reply = QMessageBox::warning(this, tr("Warning"),
                         QString("Thread number is too large! Based on the memory size, up to %1 threads can run on this computer. Are you sure you want to proceed?").arg(max_thread).arg(max_thread),
                         QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::No) {
      return true;
    }
  }
  if (ui->threadNumberEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Thread number is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->maxGapLength->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Maximum gap length is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->maxVarPTMGap->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Maximum variable PTM in a gap is empty!"),
                         QMessageBox::Yes);
    return true;
  }
  return false;
}


void TopmgWindow::updateMsg(std::string msg) {
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

void TopmgWindow::showArguments() {
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
                                        "\nshiftNumber:" + arguments_["shiftNumber"] +
                                        "\nerrorTolerance:" + arguments_["massErrorTolerance"] +
                                        "\nformErrorTolerance:" + arguments_["proteoformErrorTolerance"] +
                                        "\ncutoffSpectralType:" + arguments_["cutoffSpectralType"] +
                                        "\ncutoffSpectralValue:" + arguments_["cutoffSpectralValue"] +
                                        "\ncutoffProteoformType:" + arguments_["cutoffProteoformType"] +
                                        "\ncutoffProteoformValue:" + arguments_["cutoffProteoformValue"] +
                                        "\nallowProtMod:" + arguments_["allowProtMod"] +
                                        "\nnumOfTopPrsms:" + arguments_["numOfTopPrsms"] +
                                        "\nmaxShiftMass:" + arguments_["maxShiftMass"] +
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

void TopmgWindow::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait)
    QCoreApplication::processEvents();
}

void TopmgWindow::on_exitButton_clicked() {
  close();
}

void TopmgWindow::on_fixedModComboBox_currentIndexChanged(int index) {
  if (index == 3) {
    ui->fixedModFileEdit->setEnabled(true);
    ui->fixedModFileButton->setEnabled(true);
  } else {
    ui->fixedModFileEdit->setEnabled(false);
    ui->fixedModFileButton->setEnabled(false);
  }
}

bool TopmgWindow::nterminalerror() {
  if (ui->NONECheckBox->isChecked() || ui->NMECheckBox->isChecked() || ui->NMEACCheckBox->isChecked() || ui->MACCheckBox->isChecked()) {
    return false;
  } else {
    QMessageBox::warning(this, tr("Warning"),
                         tr("At least one N-terminal form should be selected!"),
                         QMessageBox::Yes);
    return true;
  }
}

void TopmgWindow::on_NONECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NONECheckBox->setChecked(true);
  }
}

void TopmgWindow::on_NMECheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMECheckBox->setChecked(true);
  }
}

void TopmgWindow::on_NMEACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->NMEACCheckBox->setChecked(true);
  }
}

void TopmgWindow::on_MACCheckBox_clicked(bool checked) {
  if (nterminalerror()) {
    ui->MACCheckBox->setChecked(true);
  }
}

void TopmgWindow::on_cutoffSpectralTypeComboBox_currentIndexChanged(int index) {
  if (index == 1 && !ui->decoyCheckBox->isChecked()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an FDR cutoff, the checkbox \"decoy database\" should be checked!"),
                         QMessageBox::Yes);
    ui->cutoffSpectralTypeComboBox->setCurrentIndex(0);
  }
}

void TopmgWindow::on_cutoffProteoformTypeComboBox_currentIndexChanged(int index) {
  if (index == 1 && !ui->decoyCheckBox->isChecked()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("To use an FDR cutoff, the checkbox \"decoy database\" should be checked!"),
                         QMessageBox::Yes);
    ui->cutoffProteoformTypeComboBox->setCurrentIndex(0);
  }
}

void TopmgWindow::on_decoyCheckBox_clicked(bool checked) {
  if (!checked && (ui->cutoffSpectralTypeComboBox->currentIndex() > 0 || ui->cutoffProteoformTypeComboBox->currentIndex() > 0)) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Because an FDR cutoff is selected, the checkbox \"decoy database\" cannot be unchecked."),
                         QMessageBox::Yes);
    ui->decoyCheckBox->setChecked(true);
  }
}

void TopmgWindow::closeEvent(QCloseEvent *event) {
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

bool TopmgWindow::continueToClose() {
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
