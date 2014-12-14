#-------------------------------------------------
#
# Project created by QtCreator 2014-09-03T16:42:53
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = fluids
TEMPLATE = app

INCLUDEPATH += $$_PRO_FILE_PWD_/eigen_3.2.2/ $$_PRO_FILE_PWD_/SOIL/src
LIBS += -lGLU -L$$_PRO_FILE_PWD_/SOIL/lib -lSOIL

SOURCES += main.cpp\
        mainwindow.cpp \
    glpanel.cpp \
    controller.cpp \
    simulation.cpp \
    simparameters.cpp \
    vectormath.cpp \
    camera.cpp \
    rotator.cpp

HEADERS  += mainwindow.h \
    glpanel.h \
    controller.h \
    simulation.h \
    simparameters.h \
    vectormath.h \
    camera.h \
    rotator.h \
    mat3d.h

FORMS    += mainwindow.ui

libSOIL.commands = cd $$_PRO_FILE_PWD_/SOIL/projects/makefile/ && make

QMAKE_EXTRA_TARGETS += libSOIL
PRE_TARGETDEPS += libSOIL

QMAKE_CXXFLAGS += -g -O2
