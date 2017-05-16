/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionAjuda;
    QWidget *centralWidget;
    QGridLayout *gridLayout;
    QGroupBox *groupBox_4;
    QLabel *label_5;
    QLineEdit *num_dias;
    QPushButton *calculate;
    QGroupBox *groupBox_3;
    QLineEdit *first_file;
    QLabel *label_3;
    QPushButton *first_ok;
    QLabel *label_2;
    QLineEdit *second_file;
    QPushButton *second_ok;
    QLabel *label;
    QLineEdit *third_file;
    QPushButton *third_ok;
    QLabel *label_4;
    QGroupBox *groupBox;
    QCheckBox *checkBox_ETc;
    QCheckBox *checkBox_ETr;
    QCheckBox *checkBox_Vazao;
    QCheckBox *checkBox_At;
    QCheckBox *checkBox_Pts;
    QCheckBox *checkBox_Dcr;
    QCheckBox *checkBox_Ia;
    QCheckBox *checkBox_Pe;
    QCheckBox *checkBox_vESD;
    QCheckBox *checkBox_vESS;
    QCheckBox *checkBox_vEB;
    QGroupBox *groupBox_2;
    QLabel *label_6;
    QLineEdit *lineEdit_kcr;
    QLabel *label_7;
    QLineEdit *lineEdit_kb;
    QLabel *label_8;
    QLineEdit *lineEdit_kss;
    QLineEdit *lineEdit_cs;
    QLineEdit *lineEdit_css;
    QLineEdit *lineEdit_cb;
    QLineEdit *lineEdit_coef;
    QLabel *label_9;
    QLabel *label_10;
    QLabel *label_11;
    QLabel *label_12;
    QCheckBox *checkBox_parametros;
    QLineEdit *lineEdit_maxn;
    QLabel *label_13;
    QLineEdit *lineEdit_peps;
    QLineEdit *lineEdit_ngs;
    QLabel *label_14;
    QLabel *label_15;
    QMenuBar *menuBar;
    QMenu *menuArquivo;
    QStatusBar *statusBar;
    QToolBar *mainToolBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(495, 611);
        MainWindow->setStyleSheet(QStringLiteral("background-color: rgb(160, 194, 255);"));
        actionAjuda = new QAction(MainWindow);
        actionAjuda->setObjectName(QStringLiteral("actionAjuda"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        centralWidget->setStyleSheet(QStringLiteral(""));
        gridLayout = new QGridLayout(centralWidget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        groupBox_4 = new QGroupBox(centralWidget);
        groupBox_4->setObjectName(QStringLiteral("groupBox_4"));
        label_5 = new QLabel(groupBox_4);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(10, 30, 171, 16));
        num_dias = new QLineEdit(groupBox_4);
        num_dias->setObjectName(QStringLiteral("num_dias"));
        num_dias->setGeometry(QRect(10, 50, 71, 16));
        calculate = new QPushButton(groupBox_4);
        calculate->setObjectName(QStringLiteral("calculate"));
        calculate->setGeometry(QRect(90, 50, 91, 16));
        calculate->setCursor(QCursor(Qt::PointingHandCursor));

        gridLayout->addWidget(groupBox_4, 6, 0, 1, 1);

        groupBox_3 = new QGroupBox(centralWidget);
        groupBox_3->setObjectName(QStringLiteral("groupBox_3"));
        groupBox_3->setMinimumSize(QSize(0, 374));
        groupBox_3->setMaximumSize(QSize(16777215, 374));
        first_file = new QLineEdit(groupBox_3);
        first_file->setObjectName(QStringLiteral("first_file"));
        first_file->setGeometry(QRect(10, 50, 131, 16));
        first_file->setStyleSheet(QStringLiteral(""));
        label_3 = new QLabel(groupBox_3);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(10, 30, 101, 16));
        first_ok = new QPushButton(groupBox_3);
        first_ok->setObjectName(QStringLiteral("first_ok"));
        first_ok->setGeometry(QRect(10, 70, 71, 16));
        first_ok->setCursor(QCursor(Qt::PointingHandCursor));
        label_2 = new QLabel(groupBox_3);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(10, 110, 111, 21));
        second_file = new QLineEdit(groupBox_3);
        second_file->setObjectName(QStringLiteral("second_file"));
        second_file->setGeometry(QRect(10, 130, 131, 16));
        second_ok = new QPushButton(groupBox_3);
        second_ok->setObjectName(QStringLiteral("second_ok"));
        second_ok->setGeometry(QRect(10, 150, 71, 16));
        second_ok->setCursor(QCursor(Qt::PointingHandCursor));
        label = new QLabel(groupBox_3);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(10, 190, 101, 16));
        third_file = new QLineEdit(groupBox_3);
        third_file->setObjectName(QStringLiteral("third_file"));
        third_file->setGeometry(QRect(10, 210, 131, 16));
        third_ok = new QPushButton(groupBox_3);
        third_ok->setObjectName(QStringLiteral("third_ok"));
        third_ok->setGeometry(QRect(10, 230, 71, 16));
        third_ok->setCursor(QCursor(Qt::PointingHandCursor));

        gridLayout->addWidget(groupBox_3, 1, 0, 5, 1);

        label_4 = new QLabel(centralWidget);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setCursor(QCursor(Qt::ArrowCursor));
        label_4->setFrameShadow(QFrame::Raised);
        label_4->setLineWidth(1);

        gridLayout->addWidget(label_4, 0, 0, 1, 2);

        groupBox = new QGroupBox(centralWidget);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        checkBox_ETc = new QCheckBox(groupBox);
        checkBox_ETc->setObjectName(QStringLiteral("checkBox_ETc"));
        checkBox_ETc->setGeometry(QRect(10, 30, 51, 16));
        checkBox_ETr = new QCheckBox(groupBox);
        checkBox_ETr->setObjectName(QStringLiteral("checkBox_ETr"));
        checkBox_ETr->setGeometry(QRect(10, 50, 41, 16));
        checkBox_Vazao = new QCheckBox(groupBox);
        checkBox_Vazao->setObjectName(QStringLiteral("checkBox_Vazao"));
        checkBox_Vazao->setGeometry(QRect(140, 90, 91, 16));
        checkBox_At = new QCheckBox(groupBox);
        checkBox_At->setObjectName(QStringLiteral("checkBox_At"));
        checkBox_At->setGeometry(QRect(10, 150, 101, 16));
        checkBox_Pts = new QCheckBox(groupBox);
        checkBox_Pts->setObjectName(QStringLiteral("checkBox_Pts"));
        checkBox_Pts->setGeometry(QRect(10, 90, 85, 16));
        checkBox_Dcr = new QCheckBox(groupBox);
        checkBox_Dcr->setObjectName(QStringLiteral("checkBox_Dcr"));
        checkBox_Dcr->setGeometry(QRect(10, 110, 121, 16));
        checkBox_Ia = new QCheckBox(groupBox);
        checkBox_Ia->setObjectName(QStringLiteral("checkBox_Ia"));
        checkBox_Ia->setGeometry(QRect(10, 130, 101, 16));
        checkBox_Pe = new QCheckBox(groupBox);
        checkBox_Pe->setObjectName(QStringLiteral("checkBox_Pe"));
        checkBox_Pe->setGeometry(QRect(10, 70, 101, 16));
        checkBox_vESD = new QCheckBox(groupBox);
        checkBox_vESD->setObjectName(QStringLiteral("checkBox_vESD"));
        checkBox_vESD->setGeometry(QRect(140, 30, 91, 16));
        checkBox_vESS = new QCheckBox(groupBox);
        checkBox_vESS->setObjectName(QStringLiteral("checkBox_vESS"));
        checkBox_vESS->setGeometry(QRect(140, 50, 85, 16));
        checkBox_vEB = new QCheckBox(groupBox);
        checkBox_vEB->setObjectName(QStringLiteral("checkBox_vEB"));
        checkBox_vEB->setGeometry(QRect(140, 70, 85, 16));

        gridLayout->addWidget(groupBox, 4, 1, 3, 1);

        groupBox_2 = new QGroupBox(centralWidget);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        label_6 = new QLabel(groupBox_2);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setGeometry(QRect(80, 30, 21, 16));
        lineEdit_kcr = new QLineEdit(groupBox_2);
        lineEdit_kcr->setObjectName(QStringLiteral("lineEdit_kcr"));
        lineEdit_kcr->setGeometry(QRect(10, 30, 61, 16));
        label_7 = new QLabel(groupBox_2);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setGeometry(QRect(80, 50, 21, 16));
        lineEdit_kb = new QLineEdit(groupBox_2);
        lineEdit_kb->setObjectName(QStringLiteral("lineEdit_kb"));
        lineEdit_kb->setGeometry(QRect(10, 50, 61, 16));
        label_8 = new QLabel(groupBox_2);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(80, 70, 21, 16));
        lineEdit_kss = new QLineEdit(groupBox_2);
        lineEdit_kss->setObjectName(QStringLiteral("lineEdit_kss"));
        lineEdit_kss->setGeometry(QRect(10, 70, 61, 16));
        lineEdit_cs = new QLineEdit(groupBox_2);
        lineEdit_cs->setObjectName(QStringLiteral("lineEdit_cs"));
        lineEdit_cs->setGeometry(QRect(10, 90, 61, 16));
        lineEdit_css = new QLineEdit(groupBox_2);
        lineEdit_css->setObjectName(QStringLiteral("lineEdit_css"));
        lineEdit_css->setGeometry(QRect(10, 110, 61, 16));
        lineEdit_cb = new QLineEdit(groupBox_2);
        lineEdit_cb->setObjectName(QStringLiteral("lineEdit_cb"));
        lineEdit_cb->setGeometry(QRect(10, 130, 61, 16));
        lineEdit_coef = new QLineEdit(groupBox_2);
        lineEdit_coef->setObjectName(QStringLiteral("lineEdit_coef"));
        lineEdit_coef->setGeometry(QRect(10, 150, 61, 16));
        label_9 = new QLabel(groupBox_2);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setGeometry(QRect(80, 90, 21, 16));
        label_10 = new QLabel(groupBox_2);
        label_10->setObjectName(QStringLiteral("label_10"));
        label_10->setGeometry(QRect(80, 110, 21, 16));
        label_11 = new QLabel(groupBox_2);
        label_11->setObjectName(QStringLiteral("label_11"));
        label_11->setGeometry(QRect(80, 130, 21, 16));
        label_12 = new QLabel(groupBox_2);
        label_12->setObjectName(QStringLiteral("label_12"));
        label_12->setGeometry(QRect(80, 150, 51, 20));
        checkBox_parametros = new QCheckBox(groupBox_2);
        checkBox_parametros->setObjectName(QStringLiteral("checkBox_parametros"));
        checkBox_parametros->setGeometry(QRect(10, 180, 151, 20));
        lineEdit_maxn = new QLineEdit(groupBox_2);
        lineEdit_maxn->setObjectName(QStringLiteral("lineEdit_maxn"));
        lineEdit_maxn->setGeometry(QRect(120, 30, 61, 16));
        label_13 = new QLabel(groupBox_2);
        label_13->setObjectName(QStringLiteral("label_13"));
        label_13->setGeometry(QRect(190, 30, 41, 16));
        lineEdit_peps = new QLineEdit(groupBox_2);
        lineEdit_peps->setObjectName(QStringLiteral("lineEdit_peps"));
        lineEdit_peps->setGeometry(QRect(120, 50, 61, 16));
        lineEdit_ngs = new QLineEdit(groupBox_2);
        lineEdit_ngs->setObjectName(QStringLiteral("lineEdit_ngs"));
        lineEdit_ngs->setGeometry(QRect(120, 70, 61, 16));
        label_14 = new QLabel(groupBox_2);
        label_14->setObjectName(QStringLiteral("label_14"));
        label_14->setGeometry(QRect(190, 50, 41, 16));
        label_15 = new QLabel(groupBox_2);
        label_15->setObjectName(QStringLiteral("label_15"));
        label_15->setGeometry(QRect(190, 70, 41, 16));

        gridLayout->addWidget(groupBox_2, 1, 1, 3, 1);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 495, 19));
        menuArquivo = new QMenu(menuBar);
        menuArquivo->setObjectName(QStringLiteral("menuArquivo"));
        MainWindow->setMenuBar(menuBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        MainWindow->setStatusBar(statusBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        QWidget::setTabOrder(checkBox_ETc, checkBox_ETr);
        QWidget::setTabOrder(checkBox_ETr, checkBox_Vazao);

        menuBar->addAction(menuArquivo->menuAction());
        menuArquivo->addAction(actionAjuda);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "LASH", 0));
        actionAjuda->setText(QApplication::translate("MainWindow", "&Ajuda", 0));
        actionAjuda->setShortcut(QApplication::translate("MainWindow", "Ctrl+H", 0));
        groupBox_4->setTitle(QApplication::translate("MainWindow", "Calcular resultados", 0));
        label_5->setText(QApplication::translate("MainWindow", "N\303\272mero de dias (M\303\241x. 5479)", 0));
        calculate->setText(QApplication::translate("MainWindow", "Calcular", 0));
        groupBox_3->setTitle(QApplication::translate("MainWindow", "Arquivos de entrada", 0));
        label_3->setText(QApplication::translate("MainWindow", "Primeiro Arquivo", 0));
        first_ok->setText(QApplication::translate("MainWindow", "Procurar", 0));
        label_2->setText(QApplication::translate("MainWindow", "Segundo Arquivo", 0));
        second_ok->setText(QApplication::translate("MainWindow", "Procurar", 0));
        label->setText(QApplication::translate("MainWindow", "Terceiro Arquivo", 0));
        third_ok->setText(QApplication::translate("MainWindow", "Procurar", 0));
        label_4->setText(QApplication::translate("MainWindow", "<html><head/><body><p align=\"center\"><span style=\" font-size:32pt; font-weight:600; font-style:italic; color:#070042;\">LASH</span></p></body></html>", 0));
        groupBox->setTitle(QApplication::translate("MainWindow", "Mostrar resultados para:", 0));
        checkBox_ETc->setText(QApplication::translate("MainWindow", "ETc", 0));
        checkBox_ETr->setText(QApplication::translate("MainWindow", "ETr", 0));
        checkBox_Vazao->setText(QApplication::translate("MainWindow", "Vaz\303\243o Total", 0));
        checkBox_At->setText(QApplication::translate("MainWindow", "Armaz. Total", 0));
        checkBox_Pts->setText(QApplication::translate("MainWindow", "Prec. Total", 0));
        checkBox_Dcr->setText(QApplication::translate("MainWindow", "Ascens\303\243o capilar", 0));
        checkBox_Ia->setText(QApplication::translate("MainWindow", "Abst. Inicial", 0));
        checkBox_Pe->setText(QApplication::translate("MainWindow", "Prec. Efetiva", 0));
        checkBox_vESD->setText(QApplication::translate("MainWindow", "Vaz\303\243o_ESD", 0));
        checkBox_vESS->setText(QApplication::translate("MainWindow", "Vaz\303\243o_ESS", 0));
        checkBox_vEB->setText(QApplication::translate("MainWindow", "Vaz\303\243o_EB", 0));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Par\303\242metros", 0));
        label_6->setText(QApplication::translate("MainWindow", "Kcr", 0));
        label_7->setText(QApplication::translate("MainWindow", "Kb", 0));
        label_8->setText(QApplication::translate("MainWindow", "Kss", 0));
        label_9->setText(QApplication::translate("MainWindow", "Cs", 0));
        label_10->setText(QApplication::translate("MainWindow", "Css", 0));
        label_11->setText(QApplication::translate("MainWindow", "Cb", 0));
        label_12->setText(QApplication::translate("MainWindow", "Coef. Ia", 0));
        checkBox_parametros->setText(QApplication::translate("MainWindow", "Usar valores padr\303\243o", 0));
        label_13->setText(QApplication::translate("MainWindow", "Maxn", 0));
        label_14->setText(QApplication::translate("MainWindow", "Peps", 0));
        label_15->setText(QApplication::translate("MainWindow", "Ngs", 0));
        menuArquivo->setTitle(QApplication::translate("MainWindow", "&Op\303\247\303\265es", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
