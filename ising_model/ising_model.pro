TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

unix {
    LIBS += -lblas -llapack -larmadillo
} else {
LIBS += -LC:/Armadillo/ -larmadillo
LIBS += -LC:/Armadillo/examples/lib_win64 -llapack_win64_MT -lblas_win64_MT
INCLUDEPATH += C:/Armadillo/include
}

SOURCES += main.cpp

HEADERS +=
