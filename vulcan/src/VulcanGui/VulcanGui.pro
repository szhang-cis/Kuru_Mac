#-------------------------------------------------
#
# Project created by QtCreator 2015-04-07T01:39:55
#
#-------------------------------------------------

QT       += core gui network

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = VulcanGui
TEMPLATE = app


SOURCES += main.cpp\
        vulcangui.cpp \
    instancemodel.cpp \
    socketinstance.cpp

HEADERS  += vulcangui.h \
    instancemodel.h \
    socketinstance.h

FORMS    += vulcangui.ui

RESOURCES += \
    resources.qrc
