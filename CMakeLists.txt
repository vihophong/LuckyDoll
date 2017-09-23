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
#include(${CMAKE_CURRENT_SOURCE_DIR}/aidafull.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/aidafullnew.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/aidapartial.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/aidapartialnew.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/aidasimple.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/aidasimplenew.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/aida2tree.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/aida2histobkg.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/aida2histocalib.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/belen2tree.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/merger.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/mergersimple.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/builddecay.cmake)
