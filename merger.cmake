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
project(merger)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Spectrum Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

ROOT_GENERATE_DICTIONARY(G__DataStructNew DataStructNew.h LINKDEF DataStructNewLinkDef.h)
ROOT_GENERATE_DICTIONARY(G__BELEN BELEN.h LINKDEF BELENLinkDef.h)
ROOT_GENERATE_DICTIONARY(G__Clover Clover.h LINKDEF CloverLinkDef.h)
ROOT_GENERATE_DICTIONARY(G__AIDAmerge AIDA.h LINKDEF AIDALinkDef.h)
ROOT_GENERATE_DICTIONARY(G__Beam Beam.h LINKDEF BeamLinkDef.h)
#---Create a shared library with geneated dictionary

add_library(BELEN SHARED BELEN.cpp G__BELEN.cxx) # Link2Dictionary!
add_library(Clover SHARED Clover.cpp G__Clover.cxx) # Link2Dictionary!
add_library(AIDAmerge SHARED AIDA.cpp G__AIDAmerge.cxx) # Link2Dictionary!
add_library(Beam SHARED Beam.cpp G__Beam.cxx) # Link2Dictionary!
add_library(DataStructNew SHARED DataStructNew.cpp G__DataStructNew.cxx) # Link2Dictionary!
target_link_libraries(BELEN ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(Clover ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(AIDAmerge ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(Beam ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(DataStructNew ${ROOT_LIBRARIES})   # Link2Dictionary!

#---Create  a main program using the library
add_executable(merger LuckyDollMerger.cpp  CommandLineInterface.cpp Merger.cpp CommandLineInterface.h Merger.h CommandLineInterface.h)
target_link_libraries(merger BELEN)
target_link_libraries(merger Clover)
target_link_libraries(merger AIDAmerge)
target_link_libraries(merger Beam)
target_link_libraries(merger DataStructNew)
