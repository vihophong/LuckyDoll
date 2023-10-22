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
project(mergerbeam)

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

ROOT_GENERATE_DICTIONARY(G__Beam1 Beam.h LINKDEF BeamLinkDef.h)
ROOT_GENERATE_DICTIONARY(G__BELEN1 BELEN.h LINKDEF BELENLinkDef.h)
ROOT_GENERATE_DICTIONARY(G__Clover1 Clover.h LINKDEF CloverLinkDef.h)
ROOT_GENERATE_DICTIONARY(G__AIDAmerge1 AIDA.h LINKDEF AIDALinkDef.h)
ROOT_GENERATE_DICTIONARY(G__DataStruct1 DataStruct.h LINKDEF DataStructLinkDef.h)
#---Create a shared library with geneated dictionary

add_library(BELEN1 SHARED BELEN.cpp G__BELEN1.cxx) # Link2Dictionary!
add_library(Clover1 SHARED Clover.cpp G__Clover1.cxx) # Link2Dictionary!
add_library(AIDAmerge1 SHARED AIDA.cpp G__AIDAmerge1.cxx) # Link2Dictionary!
add_library(Beam1 SHARED Beam.cpp G__Beam1.cxx) # Link2Dictionary!
add_library(DataStruct1 SHARED DataStruct.cpp G__DataStruct1.cxx) # Link2Dictionary!
target_link_libraries(BELEN1 ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(Clover1 ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(AIDAmerge1 ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(Beam1 ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(DataStruct1 ${ROOT_LIBRARIES})   # Link2Dictionary!

#---Create  a main program using the library
add_executable(mergerbeam LuckyDollMergerBeam.cpp  CommandLineInterface.cpp Merger.cpp CommandLineInterface.h Merger.h CommandLineInterface.h)
target_link_libraries(mergerbeam BELEN1)
target_link_libraries(mergerbeam Clover1)
target_link_libraries(mergerbeam Beam1)
target_link_libraries(mergerbeam AIDAmerge1)
target_link_libraries(mergerbeam DataStruct1)

