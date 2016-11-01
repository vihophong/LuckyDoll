CONFIG += c++11

TARGET = aida
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TEMPLATE = app


SOURCES += \
    #LuckyDoll.cpp \
    CloverDictionary.cpp \
    Clover.cpp \
    BELENDictionary.cpp \
    BELEN.cpp \
    AIDA.cpp \
    AIDADictionary.cpp \
    CommandLineInterface.cpp \
    AIDAUnpacker.cpp \
    BuildAIDAEvents.cpp \
    BelenReader.cpp \
    #TreeDataDictionary.cpp \
    #LuckyDoll_rawBelen.cpp \
    LuckyDollMergerSimple.cpp
    #MakeCalibHisto.cpp

HEADERS += \
    TreeData.h \
    CloverDictionary.h \
    Clover.h \
    BELENdefs.h \
    BELENDictionary.h \
    BELEN.h \
    AIDA.h \
    AIDAdefs.h \
    AIDALinkDef.h \
    AIDADictionary.h \
    CommandLineInterface.h \
    rawaida.h \
    AIDAUnpacker.h \
    BuildAIDAEvents.h \
    BelenReader.h \
    TreeDataLinkDef.h \
    #TreeDataDictionary.h

#---------------------ROOT include----------------------------
LIBS += $$system(root-config --cflags --glibs)  -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl

INCLUDEPATH += "$(ROOTSYS)/include"
DEPENDPATH += "$(ROOTSYS)/include"

#---------------------ANAROOT include----------------------------
LIBS +=-L$(TARTSYS)/lib -lanaeurica -lananadeko -lanacore -lanabrips -lXMLParser
INCLUDEPATH +=$(TARTSYS)/include

