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

#include "gui/util/gui_check.hpp"

#include <QMessageBox>
#include <QObject>

namespace toppic {

namespace gui_check {

bool invalidEdit(QWidget* parent, const QLineEdit* edit, const QString& name) {
  QString msg;
  if (edit->text().trimmed().isEmpty()) {
    msg = QObject::tr("%1 is empty!").arg(name);
  } else if (!edit->hasAcceptableInput()) {
    msg = QObject::tr("%1 %2 is not a valid value!").arg(name, edit->text());
  } else {
    return false;
  }
  QMessageBox::warning(parent, QObject::tr("Warning"), msg, QMessageBox::Yes);
  return true;
}

bool pathWithSpace(QWidget* parent, const QString& path, const QString& name) {
  if (!path.contains(' ')) {
    return false;
  }
  QMessageBox::warning(parent, QObject::tr("Warning"),
                       QObject::tr("The %1 %2 contains a space! Rename it or "
                                   "move it to a folder without spaces.")
                           .arg(name, path),
                       QMessageBox::Yes);
  return true;
}

}  // namespace gui_check

}  // namespace toppic
