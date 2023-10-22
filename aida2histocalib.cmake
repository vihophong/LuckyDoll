#####################################################################################################################

# Jorge Agramunt Ros    @ IFIC(Valencia,Spain)  jorge.agramunt@ific.uv.es
# Alvaro Tolosa Delgado @ IFIC(Valencia,Spain)  alvaro.tolosa@ific.uv.es
# Copyright (c) 2016 Jorge Agramunt & Alvaro Tolosa. All rights reserved.
#####################################################################################################################

# Example adapted from https://root.cern.ch/faq/can-i-integrate-root-my-cmake-build

#####################################################################################################################

# For compiling:
# mkdir build
# cd build
# cmake ..
# make

#####################################################################################################################


# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(aidacalib)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

ROOT_GENERATE_DICTIONARY(G__AIDAcal AIDA.h LINKDEF AIDALinkDef.h)

#---Create a shared library with geneated dictionary
add_library(AIDAcal SHARED AIDA.cpp G__AIDAcal.cxx) # Link2Dictionary!
target_link_libraries(AIDAcal ${ROOT_LIBRARIES})   # Link2Dictionary!

#boost lib
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
#minimum required version of boost is 1.42
find_package(Boost 1.42.0 COMPONENTS system iostreams)
#include_directories(${Boost_INCLUDE_DIRS})

#---Create  a main program using the library
add_executable(aida2histcalib MakeCalibHisto.cpp AIDAUnpackerGz.cpp BuildAIDAEventsNew.cpp  CommandLineInterface.cpp AIDAUnpackerGz.h BuildAIDAEventsNew.h CommandLineInterface.h rawaida.h)
target_link_libraries(aida2histcalib AIDAcal)
target_link_libraries(aida2histcalib ${Boost_LIBRARIES})
