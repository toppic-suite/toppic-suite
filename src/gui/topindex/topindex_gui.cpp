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

#include <QLocale>
#include <QApplication>
#include <QGuiApplication>
#include <QScreen>
#include <QFontDatabase>
#include <QDesktopWidget>

#include "gui/topindex/topindexdialog.hpp"

int main(int argc, char *argv[]) {
  // make sure we are using the c locale
  QLocale::setDefault(QLocale::c());
  QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

  QApplication app(argc, argv);

  QScreen *screen = QGuiApplication::primaryScreen();
  QRect deskRect = screen->availableGeometry();

  TopIndexDialog td;
  td.show();
  td.move((deskRect.width() - td.width()) / 2, (deskRect.height() - td.height()) / 2);

  return app.exec();
}
