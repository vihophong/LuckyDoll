CONFIG += c++11

TARGET = aida
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TEMPLATE = app


SOURCES += \
    LuckyDoll.cpp \
    AIDA.cpp \
    AIDADictionary.cpp \
    CommandLineInterface.cpp \
    AIDAUnpacker.cpp \
    BuildAIDAEvents.cpp \
    #MakeBackgroundHisto.cpp

HEADERS += \
    AIDA.h \
    AIDAdefs.h \
    AIDALinkDef.h \
    AIDADictionary.h \
    CommandLineInterface.h \
    rawaida.h \
    AIDAUnpacker.h \
    BuildAIDAEvents.h

#---------------------ROOT include----------------------------
LIBS += $$system(root-config --cflags --glibs)  -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl

INCLUDEPATH += "$(ROOTSYS)/include"
DEPENDPATH += "$(ROOTSYS)/include"

#---------------------ANAROOT include----------------------------
LIBS +=-L$(TARTSYS)/lib -lanaeurica -lananadeko -lanacore -lanabrips -lXMLParser
INCLUDEPATH +=$(TARTSYS)/include

#---lib local
LIBS +=-L/home/phong/lib -lAIDA

