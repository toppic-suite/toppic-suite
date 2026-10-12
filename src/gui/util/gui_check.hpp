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

#ifndef TOPPIC_GUI_UTIL_GUI_CHECK_HPP_
#define TOPPIC_GUI_UTIL_GUI_CHECK_HPP_

#include <QLineEdit>
#include <QString>
#include <QWidget>

namespace toppic {

// Checks of the GUI fields before a tool is started; each shows a warning
// and returns true when the check fails.
namespace gui_check {

// The field is empty or its text is not accepted by its validator.
bool invalidEdit(QWidget* parent, const QLineEdit* edit, const QString& name);

// The path contains a space: the GUIs pass the command line to QProcess split
// at spaces, so such a path would be cut in two.
bool pathWithSpace(QWidget* parent, const QString& path, const QString& name);

}  // namespace gui_check

}  // namespace toppic

#endif
