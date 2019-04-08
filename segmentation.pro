#-------------------------------------------------
#
# Project created by QtCreator
#
#-------------------------------------------------

QT       += core gui
QT       += opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = segmentation
TEMPLATE = app

DEFINES += QT_DEPRECATED_WARNINGS

LIBS += -lglut -lGLU
LIBS += -L$$PWD/../OpenMesh/liblinux/ -lOpenMeshCore

INCLUDEPATH += $$PWD/../OpenMesh/inc/
DEPENDPATH += $$PWD/../OpenMesh/inc/
DEPENDPATH += $$PWD/../OpenMesh/liblinux/


SOURCES += \
        main.cpp \
        mainwindow.cpp \
    meshviewerwidget.cpp

HEADERS += \
        mainwindow.h \
    meshviewerwidget.h

FORMS += \
        mainwindow.ui


