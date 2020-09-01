QT += opengl
TEMPLATE = app
TARGET = a.out
INCLUDEPATH += .

DEFINES += QT_DEPRECATED_WARNINGS

HEADERS += glwidget.h init.h solve.h approximation.h
SOURCES += main.cpp glwidget.cpp init.cpp solve.cpp approximation.cpp 
