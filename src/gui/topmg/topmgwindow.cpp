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


#include <QFileDialog>
#include <QElapsedTimer>
#include <QMessageBox>
#include <QCloseEvent>
#include <QCoreApplication>
#include <QToolTip>
#include <QDesktopServices>

#include "topmgwindow.h"
#include "ui_topmgwindow.h"
#include "threadtopmg.h"

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

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
      font.setFamily(QStringLiteral("Courier New"));
      fontTable.setFamily(QStringLiteral("Courier New"));
#else
      font.setFamily(QStringLiteral("Monospace"));
      fontTable.setFamily(QStringLiteral("Monospace"));
#endif
      ui->outputTextBrowser->setFont(font);
      thread_ = new threadtopmg(this);
      showInfo = "";
      setToolTip("");
      setToolTipDuration(100);

      fontTable.setPointSize(9);
      ui->listWidget->setFont(fontTable);

      on_defaultButton_clicked();
    }

topmgWindow::~topmgWindow() {
  delete ui;
}

void topmgWindow::initArguments() {
  arguments_["oriDatabaseFileName"]="";
  arguments_["databaseFileName"] = "";
  arguments_["databaseBlockSize"] = "1000000";
  arguments_["spectrumFileName"] = "";
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
  arguments_["useGf"] = "false";
  arguments_["executiveDir"] = ".";
  arguments_["resourceDir"] = "";
  arguments_["keepTempFiles"] = "false";
  arguments_["fullBinaryPath"] = "false";
  arguments_["groupSpectrumNumber"] = "1";
  arguments_["filteringResultNumber"] = "20";
  arguments_["residueModFileName"] = "";
  arguments_["threadNumber"] = "1";
  arguments_["useFeatureFile"] = "true";
  arguments_["skipList"] = "";
  arguments_["proteo_graph_dis"] = "40";
  arguments_["useASFDiag"] = "false";
  arguments_["varPtmNumber"] = "10";
  arguments_["varPtmNumInGap"] = "5";
}

void topmgWindow::on_clearButton_clicked() {
  ui->databaseFileEdit->clear();
  ui->fixedModFileEdit->clear();
  ui->errorToleranceEdit->setText("15");
  ui->maxModEdit->clear();
  ui->cutoffSpectralValueEdit->clear();
  ui->cutoffProteoformValueEdit->clear();
  ui->modFileEdit->clear();
  ui->threadNumberEdit->clear();
  ui->skipListEdit->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum file.");
  ui->fixedModComboBox->setCurrentIndex(0);
  on_fixedModComboBox_currentIndexChanged(0);
  ui->activationComboBox->setCurrentIndex(0);
  ui->cutoffSpectralTypeComboBox->setCurrentIndex(0);
  ui->cutoffProteoformTypeComboBox->setCurrentIndex(0);
  ui->numModComboBox->setCurrentIndex(4);
  ui->numUnknownShiftComboBox->setCurrentIndex(0);
  ui->NONECheckBox->setChecked(false);
  ui->NMECheckBox->setChecked(false);
  ui->NMEACCheckBox->setChecked(false);
  ui->MACCheckBox->setChecked(false);
  ui->decoyCheckBox->setChecked(false);
  ui->topfdFeatureCheckBox->setChecked(false);
}

void topmgWindow::on_defaultButton_clicked() {
  ui->fixedModFileEdit->clear();
  ui->errorToleranceEdit->setText("15");
  ui->maxModEdit->setText("500");
  ui->cutoffSpectralValueEdit->setText("0.01");
  ui->cutoffProteoformValueEdit->setText("0.01");
  ui->threadNumberEdit->setText("1");
  ui->skipListEdit->clear();
  ui->outputTextBrowser->setText("Click the Start button to process the spectrum file.");
  ui->fixedModComboBox->setCurrentIndex(0);
  on_fixedModComboBox_currentIndexChanged(0);
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
}

void topmgWindow::updatedir(QString s) {
  if (!s.isEmpty()) {
    lastDir_ = s;
  }
}

// single file
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
      "Fixed modification files(*.txt);;All files(*.*)");
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
      "Modification files(*.txt);;All files(*.*)");
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
  std::string lastinfo = "";
  std::string nowinfo = "";
  std::string info;
  while (true) {
    // Here is the infomation been shown in the infoBox.
    nowinfo = buffer.str();
    info = nowinfo.substr(lastinfo.length());
    lastinfo = nowinfo;
    if (info != "") {
      updateMsg(info);
    }
    if (thread_->isFinished()) {
      break;
    }
    sleep(10);
  }
  unlockDialog();

  showInfo = "";
  thread_->exit();
  std::cout.rdbuf(oldbuf);
}

void topmgWindow::on_outputButton_clicked() {
  fs::path full_path(arguments_["spectrumFileName"].c_str());
  QString outPath = full_path.remove_filename().string().c_str();
  QDesktopServices::openUrl(QUrl(outPath, QUrl::TolerantMode));
}

std::map<std::string, std::string> topmgWindow::getArguments() {


  QString path = QCoreApplication::applicationFilePath();
  std::string exe_dir = prot::file_util::getExecutiveDir(path.toStdString());
  arguments_["executiveDir"] = exe_dir;
  arguments_["resourceDir"] = arguments_["executiveDir"] + prot::file_util::getFileSeparator() + prot::file_util::getResourceDirName();
  arguments_["oriDatabaseFileName"] = ui->databaseFileEdit->text().toStdString();
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
  arguments_["fullBinaryPath"] = "false";  // default
  arguments_["filteringResultNumber"] = "20";  // default
  arguments_["useGf"] = "false";  // default
  arguments_["groupSpectrumNumber"] = "1";  // default
  arguments_["skipList"] = ui->skipListEdit->text().toStdString();
  arguments_["proteo_graph_dis"] = "40";  // default
  arguments_["useASFDiag"] = "false";  // default
  arguments_["varPtmNumber"] = ui->numModComboBox->currentText().toStdString();
  arguments_["ptmNumber"] = ui->numUnknownShiftComboBox->currentText().toStdString();
  arguments_["varPtmNumInGap"] = "5";  // default
  arguments_["residueModFileName"] = ui->modFileEdit->text().toStdString();
  arguments_["threadNumber"] = ui->threadNumberEdit->text().toStdString();
  if (ui->topfdFeatureCheckBox->isChecked()) {
    arguments_["useFeatureFile"] = "false";
  } else {
    arguments_["useFeatureFile"] = "true";
  }
  // showArguments();
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
      "Spectrum files(*.msalign)");
  for (int i = 0; i < spfiles.size(); i++) {
    QString spfile = spfiles.at(i);
    updatedir(spfile);
    if (ableToAdd(spfile)) {
      ui->listWidget->addItem(new QListWidgetItem(spfile));
    }

  }
};

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
};

void topmgWindow::lockDialog() {
  ui->databaseFileButton->setEnabled(false);
  ui->modFileButton->setEnabled(false);
  ui->databaseFileEdit->setEnabled(false);
  ui->fixedModFileEdit->setEnabled(false);
  ui->errorToleranceEdit->setEnabled(false);
  ui->maxModEdit->setEnabled(false);
  ui->cutoffSpectralValueEdit->setEnabled(false);
  ui->cutoffProteoformValueEdit->setEnabled(false);
  ui->modFileEdit->setEnabled(false);
  ui->threadNumberEdit->setEnabled(false);
  ui->skipListEdit->setEnabled(false);
  ui->fixedModComboBox->setEnabled(false);
  on_fixedModComboBox_currentIndexChanged(0);
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
  if (ui->modFileEdit->text().isEmpty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("Please select a modification file for variable PTMs!"),
                         QMessageBox::Yes);
    return true;
  }
  return false;
}

void topmgWindow::updateMsg(std::string msg) {
  QString info = msg.c_str();
  int lastloc = info.lastIndexOf("\r", info.size() - 2);
  if (lastloc > 0) {
    info = info.right(info.size() - lastloc - 1);
  }
  if (info.at(0) == '\n') {
    info = info.right(info.size() - 1);
  }
  ui->outputTextBrowser->setText(showInfo + info);
  if (msg.at(msg.size() - 1) != '\r' && info != "\n") {
    showInfo = ui->outputTextBrowser->toPlainText();
  }
  QTextCursor cursor = ui->outputTextBrowser->textCursor();
  cursor.movePosition(QTextCursor::End);
  ui->outputTextBrowser->setTextCursor(cursor);
}

void topmgWindow::showArguments() {
  QMessageBox::warning(0, "Arguments", ("executiveDir:" + arguments_["executiveDir"] +
                                        "\nresourceDir:" + arguments_["resourceDir"] +
                                        "\noriDatabaseFileName:" + arguments_["oriDatabaseFileName"] +
                                        "\ndatabaseFileName:" + arguments_["databaseFileName"] +
                                        "\ndatabaseBlockSize:" + arguments_["databaseBlockSize"] +
                                        "\nspectrumFileName:" + arguments_["combinedOutputName"] +
                                        "\ncombinedOutputName:" + arguments_["spectrumFileName"] +
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
                                        "\nuseGf:" + arguments_["useGf"] +
                                        "\nkeepTempFiles:" + arguments_["keepTempFiles"] +
                                        "\nfullBinaryPath:" + arguments_["fullBinaryPath"] +
                                        "\nlocal_threshold:" + arguments_["local_threshold"] +
                                        "\ngroupSpectrumNumber:" + arguments_["groupSpectrumNumber"] +
                                        "\nfilteringResultNumber:" + arguments_["filteringResultNumber"] +
                                        "\nresidueModFileName:" + arguments_["residueModFileName"] +
                                        "\nthreadNumber:" + arguments_["threadNumber"] +
                                        "\nuseFeatureFile:" + arguments_["useFeatureFile"] +
                                        "\nskipList:" + arguments_["skipList"]).c_str(), QMessageBox::Yes);
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
