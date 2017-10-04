/********************************************************************************
** Form generated from reading UI file 'topfddialog.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TOPFDDIALOG_H
#define UI_TOPFDDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_TopFDDialog
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QLineEdit *fileEdit;
    QPushButton *fileButton;
    QWidget *horizontalLayoutWidget_2;
    QHBoxLayout *hboxLayout;
    QLabel *label_2;
    QLineEdit *maxChargeEdit;
    QWidget *horizontalLayoutWidget_3;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_3;
    QLineEdit *maxMassEdit;
    QWidget *horizontalLayoutWidget_4;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_4;
    QLineEdit *mzErrorEdit;
    QWidget *horizontalLayoutWidget_5;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_5;
    QLineEdit *snRatioEdit;
    QWidget *horizontalLayoutWidget_6;
    QHBoxLayout *horizontalLayout_6;
    QCheckBox *missLevelOneCheckBox;
    QProgressBar *progressBar;
    QWidget *horizontalLayoutWidget_7;
    QHBoxLayout *horizontalLayout_7;
    QPushButton *exitButton;
    QPushButton *clearButton;
    QPushButton *defaultButton;
    QPushButton *startButton;
    QGroupBox *msgBox;
    QLabel *infolabel;

    void setupUi(QDialog *TopFDDialog)
    {
        if (TopFDDialog->objectName().isEmpty())
            TopFDDialog->setObjectName(QStringLiteral("TopFDDialog"));
        TopFDDialog->resize(578, 353);
        horizontalLayoutWidget = new QWidget(TopFDDialog);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(31, 20, 511, 31));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(horizontalLayoutWidget);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout->addWidget(label);

        fileEdit = new QLineEdit(horizontalLayoutWidget);
        fileEdit->setObjectName(QStringLiteral("fileEdit"));
        fileEdit->setReadOnly(true);

        horizontalLayout->addWidget(fileEdit);

        fileButton = new QPushButton(horizontalLayoutWidget);
        fileButton->setObjectName(QStringLiteral("fileButton"));

        horizontalLayout->addWidget(fileButton);

        horizontalLayoutWidget_2 = new QWidget(TopFDDialog);
        horizontalLayoutWidget_2->setObjectName(QStringLiteral("horizontalLayoutWidget_2"));
        horizontalLayoutWidget_2->setGeometry(QRect(30, 60, 203, 31));
        hboxLayout = new QHBoxLayout(horizontalLayoutWidget_2);
        hboxLayout->setObjectName(QStringLiteral("hboxLayout"));
        hboxLayout->setSizeConstraint(QLayout::SetDefaultConstraint);
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(horizontalLayoutWidget_2);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setMinimumSize(QSize(120, 0));
        label_2->setLayoutDirection(Qt::LeftToRight);

        hboxLayout->addWidget(label_2);

        maxChargeEdit = new QLineEdit(horizontalLayoutWidget_2);
        maxChargeEdit->setObjectName(QStringLiteral("maxChargeEdit"));
        maxChargeEdit->setEnabled(true);
        maxChargeEdit->setMinimumSize(QSize(75, 0));
        maxChargeEdit->setMaximumSize(QSize(70, 16777215));

        hboxLayout->addWidget(maxChargeEdit);

        horizontalLayoutWidget_3 = new QWidget(TopFDDialog);
        horizontalLayoutWidget_3->setObjectName(QStringLiteral("horizontalLayoutWidget_3"));
        horizontalLayoutWidget_3->setGeometry(QRect(30, 100, 203, 31));
        horizontalLayout_3 = new QHBoxLayout(horizontalLayoutWidget_3);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        label_3 = new QLabel(horizontalLayoutWidget_3);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setMinimumSize(QSize(120, 0));

        horizontalLayout_3->addWidget(label_3);

        maxMassEdit = new QLineEdit(horizontalLayoutWidget_3);
        maxMassEdit->setObjectName(QStringLiteral("maxMassEdit"));
        maxMassEdit->setMinimumSize(QSize(75, 0));
        maxMassEdit->setMaximumSize(QSize(70, 16777215));

        horizontalLayout_3->addWidget(maxMassEdit);

        horizontalLayoutWidget_4 = new QWidget(TopFDDialog);
        horizontalLayoutWidget_4->setObjectName(QStringLiteral("horizontalLayoutWidget_4"));
        horizontalLayoutWidget_4->setGeometry(QRect(350, 60, 193, 31));
        horizontalLayout_4 = new QHBoxLayout(horizontalLayoutWidget_4);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, 0, 0, 0);
        label_4 = new QLabel(horizontalLayoutWidget_4);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setMinimumSize(QSize(110, 0));

        horizontalLayout_4->addWidget(label_4);

        mzErrorEdit = new QLineEdit(horizontalLayoutWidget_4);
        mzErrorEdit->setObjectName(QStringLiteral("mzErrorEdit"));
        mzErrorEdit->setMinimumSize(QSize(75, 0));
        mzErrorEdit->setMaximumSize(QSize(70, 16777215));

        horizontalLayout_4->addWidget(mzErrorEdit);

        horizontalLayoutWidget_5 = new QWidget(TopFDDialog);
        horizontalLayoutWidget_5->setObjectName(QStringLiteral("horizontalLayoutWidget_5"));
        horizontalLayoutWidget_5->setGeometry(QRect(350, 100, 193, 31));
        horizontalLayout_5 = new QHBoxLayout(horizontalLayoutWidget_5);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(0, 0, 0, 0);
        label_5 = new QLabel(horizontalLayoutWidget_5);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setMinimumSize(QSize(110, 0));

        horizontalLayout_5->addWidget(label_5);

        snRatioEdit = new QLineEdit(horizontalLayoutWidget_5);
        snRatioEdit->setObjectName(QStringLiteral("snRatioEdit"));
        snRatioEdit->setMinimumSize(QSize(75, 0));
        snRatioEdit->setMaximumSize(QSize(70, 16777215));

        horizontalLayout_5->addWidget(snRatioEdit);

        horizontalLayoutWidget_6 = new QWidget(TopFDDialog);
        horizontalLayoutWidget_6->setObjectName(QStringLiteral("horizontalLayoutWidget_6"));
        horizontalLayoutWidget_6->setGeometry(QRect(30, 140, 160, 31));
        horizontalLayout_6 = new QHBoxLayout(horizontalLayoutWidget_6);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        horizontalLayout_6->setContentsMargins(0, 0, 0, 0);
        missLevelOneCheckBox = new QCheckBox(horizontalLayoutWidget_6);
        missLevelOneCheckBox->setObjectName(QStringLiteral("missLevelOneCheckBox"));
        missLevelOneCheckBox->setAcceptDrops(false);
        missLevelOneCheckBox->setCheckable(true);
        missLevelOneCheckBox->setChecked(false);

        horizontalLayout_6->addWidget(missLevelOneCheckBox);

        progressBar = new QProgressBar(TopFDDialog);
        progressBar->setObjectName(QStringLiteral("progressBar"));
        progressBar->setGeometry(QRect(30, 270, 511, 23));
        progressBar->setAutoFillBackground(false);
        progressBar->setValue(0);
        progressBar->setTextVisible(true);
        horizontalLayoutWidget_7 = new QWidget(TopFDDialog);
        horizontalLayoutWidget_7->setObjectName(QStringLiteral("horizontalLayoutWidget_7"));
        horizontalLayoutWidget_7->setGeometry(QRect(30, 300, 511, 41));
        horizontalLayout_7 = new QHBoxLayout(horizontalLayoutWidget_7);
        horizontalLayout_7->setSpacing(30);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        horizontalLayout_7->setContentsMargins(0, 0, 0, 0);
        exitButton = new QPushButton(horizontalLayoutWidget_7);
        exitButton->setObjectName(QStringLiteral("exitButton"));

        horizontalLayout_7->addWidget(exitButton);

        clearButton = new QPushButton(horizontalLayoutWidget_7);
        clearButton->setObjectName(QStringLiteral("clearButton"));

        horizontalLayout_7->addWidget(clearButton);

        defaultButton = new QPushButton(horizontalLayoutWidget_7);
        defaultButton->setObjectName(QStringLiteral("defaultButton"));

        horizontalLayout_7->addWidget(defaultButton);

        startButton = new QPushButton(horizontalLayoutWidget_7);
        startButton->setObjectName(QStringLiteral("startButton"));
        startButton->setEnabled(true);
        startButton->setAcceptDrops(false);
        startButton->setCheckable(false);
        startButton->setAutoDefault(true);
        startButton->setFlat(false);

        horizontalLayout_7->addWidget(startButton);

        msgBox = new QGroupBox(TopFDDialog);
        msgBox->setObjectName(QStringLiteral("msgBox"));
        msgBox->setGeometry(QRect(30, 180, 511, 71));
        msgBox->setStyleSheet(QLatin1String("QGroupBox{\n"
"border-width:1px;\n"
"border-style:solid;\n"
"border-radius: 5px;\n"
"border-color:gray;\n"
"margin-top:1.5ex;\n"
"}\n"
"QGroupBox::title{\n"
"subcontrol-origin:margin;\n"
"subcontrol-position:top left;\n"
"left:10px;\n"
"margin-left:0px;\n"
"padding:0 1px;\n"
"}"));
        infolabel = new QLabel(msgBox);
        infolabel->setObjectName(QStringLiteral("infolabel"));
        infolabel->setEnabled(true);
        infolabel->setGeometry(QRect(10, 30, 511, 31));
        infolabel->setFrameShape(QFrame::NoFrame);
        infolabel->setTextFormat(Qt::AutoText);
        infolabel->setWordWrap(true);
        infolabel->setOpenExternalLinks(false);

        retranslateUi(TopFDDialog);
        QObject::connect(exitButton, SIGNAL(clicked()), TopFDDialog, SLOT(close()));

        startButton->setDefault(false);


        QMetaObject::connectSlotsByName(TopFDDialog);
    } // setupUi

    void retranslateUi(QDialog *TopFDDialog)
    {
        TopFDDialog->setWindowTitle(QApplication::translate("TopFDDialog", "TopFD", 0));
        label->setText(QApplication::translate("TopFDDialog", "Spectra file : ", 0));
        fileButton->setText(QApplication::translate("TopFDDialog", "File", 0));
        label_2->setText(QApplication::translate("TopFDDialog", "Max charge (Da) : ", 0));
        maxChargeEdit->setText(QApplication::translate("TopFDDialog", "30", 0));
        label_3->setText(QApplication::translate("TopFDDialog", "Max mass (Da) : ", 0));
        maxMassEdit->setText(QApplication::translate("TopFDDialog", "100000", 0));
        label_4->setText(QApplication::translate("TopFDDialog", "Mz error (Da) : ", 0));
        mzErrorEdit->setText(QApplication::translate("TopFDDialog", "0.02", 0));
        label_5->setText(QApplication::translate("TopFDDialog", "S/N ratio : ", 0));
        snRatioEdit->setText(QApplication::translate("TopFDDialog", "1.0", 0));
        missLevelOneCheckBox->setText(QApplication::translate("TopFDDialog", "Missing level one", 0));
        exitButton->setText(QApplication::translate("TopFDDialog", "Exit", 0));
        clearButton->setText(QApplication::translate("TopFDDialog", "Clear", 0));
        defaultButton->setText(QApplication::translate("TopFDDialog", "Default", 0));
        startButton->setText(QApplication::translate("TopFDDialog", "Start", 0));
        msgBox->setTitle(QApplication::translate("TopFDDialog", "Message box", 0));
        infolabel->setText(QApplication::translate("TopFDDialog", "Press start to process the speatra.", 0));
    } // retranslateUi

};

namespace Ui {
    class TopFDDialog: public Ui_TopFDDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TOPFDDIALOG_H
