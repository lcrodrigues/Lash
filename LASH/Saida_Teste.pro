#-------------------------------------------------
#
# Project created by QtCreator 2016-08-18T04:20:57
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = lash_ui
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    varclasses.cpp \
    secondwindow.cpp \
    helpdialog.cpp

HEADERS  += mainwindow.h \
    varclasses.h \
    secondwindow.h \
    helpdialog.h

FORMS    += mainwindow.ui \
    secondwindow.ui \
    helpdialog.ui
