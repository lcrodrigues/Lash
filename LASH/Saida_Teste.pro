#-------------------------------------------------
#
# Project created by QtCreator 2016-08-18T04:20:57
#
#-------------------------------------------------

QT       += core gui
QT       += charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = lash_ui
TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

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
