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

#include "gui/top_converter/topconverterdialog.hpp"

#include <QCloseEvent>
#include <QDesktopServices>
#include <QDoubleValidator>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QFontDatabase>
#include <QMessageBox>
#include <QScrollBar>
#include <QUrl>
#include <string>

#include "common/util/file_util.hpp"
#include "common/util/version.hpp"
#include "gui/top_converter/ui_topconverterdialog.h"
#include "gui/util/command.hpp"
#include "gui/util/gui_message.hpp"
#include "sql/ms_sql/ms_sql_writer.hpp"

namespace {

const char* kIdleMsg = "Click the Start button to convert the spectrum file.";

}  // namespace

TopConverterDialog::TopConverterDialog(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::TopConverterDialog) {
  ui->setupUi(this);
  std::string title = "TopConverter v." + toppic::Version::getVersion();
  this->setWindowTitle(QString::fromStdString(title));
  lastDir_ = ".";

  ui->mzSizeEdit->setValidator(new QDoubleValidator(0.0, 1e6, 6, this));
  ui->rtDividerEdit->setValidator(new QDoubleValidator(0.0, 1e6, 6, this));

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

  TopConverterDialog::on_defaultButton_clicked();
}

TopConverterDialog::~TopConverterDialog() {
  if (process_.state() != QProcess::NotRunning) {
    process_.kill();
  }
  delete ui;
}

void TopConverterDialog::closeEvent(QCloseEvent* event) {
  if (process_.state() != QProcess::NotRunning) {
    if (!continueToClose()) {
      event->ignore();
      return;
    }
    process_.kill();
  }
  event->accept();
}

void TopConverterDialog::on_browseButton_clicked() {
  QString file = QFileDialog::getOpenFileName(
      this, "Select a spectrum file", lastDir_,
      "Spectrum files (*.mzML *.mzXML)");
  if (file.isEmpty()) {
    return;
  }
  if (file.toStdString().length() > 200) {
    QMessageBox::warning(this, tr("Warning"), tr("The file path is too long!"),
                         QMessageBox::Yes);
    return;
  }
  ui->spectrumFileEdit->setText(file);
  lastDir_ = QString::fromStdString(
      toppic::file_util::directory(file.toStdString()));
}

void TopConverterDialog::on_clearButton_clicked() {
  ui->spectrumFileEdit->clear();
  ui->outputTextBrowser->setText(kIdleMsg);
  ui->outputButton->setEnabled(false);
}

void TopConverterDialog::on_defaultButton_clicked() {
  ui->mzSizeEdit->setText(
      QString::number(toppic::MsSqlWriter::DEFAULT_MZ_SIZE));
  ui->rtDividerEdit->setText(
      QString::number(toppic::MsSqlWriter::DEFAULT_RT_DIVIDER));
  ui->outputTextBrowser->setText(kIdleMsg);
}

std::string TopConverterDialog::getSpectrumFileName() {
  return ui->spectrumFileEdit->text().trimmed().toStdString();
}

void TopConverterDialog::on_startButton_clicked() {
  if (checkError()) {
    return;
  }
  lockDialog();

  QString path = QCoreApplication::applicationFilePath();
  std::string exe_dir = toppic::file_util::getExecutiveDir(path.toStdString());
  if (toppic::file_util::checkSpace(exe_dir)) {
    ui->outputTextBrowser->setText(
        "Current directory " + QString::fromStdString(exe_dir) +
        " contains space and will cause errors in the program!");
  }
  std::string cmd = toppic::command::geneTopConverterCommand(
      exe_dir, getSpectrumFileName(),
      ui->mzSizeEdit->text().toStdString(),
      ui->rtDividerEdit->text().toStdString());
  QString q_cmd = QString::fromStdString(cmd).trimmed();
  QStringList cmd_list = q_cmd.split(" ");
  QString prog = cmd_list[0];
  cmd_list.removeFirst();

  process_.start(prog, cmd_list);
  process_.waitForStarted();

  toppic::GuiMessage guiMsg;
  bool finish = false;
  while (!finish) {
    if (process_.state() == QProcess::NotRunning) {
      finish = true;
    }
    bool ready = process_.waitForReadyRead(100);
    if (ready || finish) {
      QByteArray byteArray = process_.readAllStandardOutput();
      std::string msg = guiMsg.getMsg(QString(byteArray).toStdString());
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

void TopConverterDialog::on_exitButton_clicked() { close(); }

bool TopConverterDialog::continueToClose() {
  return QMessageBox::question(
             this, tr("Quit"),
             tr("TopConverter is still running. Are you sure you want to "
                "quit?"),
             QMessageBox::Yes | QMessageBox::No,
             QMessageBox::No) == QMessageBox::Yes;
}

void TopConverterDialog::on_outputButton_clicked() {
  std::string dir = toppic::file_util::directory(getSpectrumFileName());
  QDesktopServices::openUrl(
      QUrl::fromLocalFile(QString::fromStdString(dir)));
}

void TopConverterDialog::lockDialog() {
  ui->browseButton->setEnabled(false);
  ui->spectrumFileEdit->setEnabled(false);
  ui->mzSizeEdit->setEnabled(false);
  ui->rtDividerEdit->setEnabled(false);
  ui->clearButton->setEnabled(false);
  ui->defaultButton->setEnabled(false);
  ui->startButton->setEnabled(false);
  ui->outputButton->setEnabled(false);
}

void TopConverterDialog::unlockDialog() {
  ui->browseButton->setEnabled(true);
  ui->spectrumFileEdit->setEnabled(true);
  ui->mzSizeEdit->setEnabled(true);
  ui->rtDividerEdit->setEnabled(true);
  ui->clearButton->setEnabled(true);
  ui->defaultButton->setEnabled(true);
  ui->startButton->setEnabled(true);
  ui->outputButton->setEnabled(true);
  ui->outputButton->setDefault(true);
}

bool TopConverterDialog::checkError() {
  if (getSpectrumFileName().empty()) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("The spectrum file is not selected!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->mzSizeEdit->text().toDouble() <= 0.0) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("The m/z block width must be positive!"),
                         QMessageBox::Yes);
    return true;
  }
  if (ui->rtDividerEdit->text().toDouble() <= 0.0) {
    QMessageBox::warning(this, tr("Warning"),
                         tr("The retention time divider must be positive!"),
                         QMessageBox::Yes);
    return true;
  }
  return false;
}

void TopConverterDialog::updateMsg(std::string msg) {
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

void TopConverterDialog::sleep(int wait) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < wait) {
    QCoreApplication::processEvents();
  }
}
