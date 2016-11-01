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
project(AllProjects)
include(${CMAKE_CURRENT_SOURCE_DIR}/aidafull.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/aidasimple.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/aida2tree.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/aida2histobkg.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/aida2histocalib.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/belen2tree.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/mergersimple.cmake)
