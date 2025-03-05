cmake_minimum_required(VERSION 3.27 FATAL_ERROR)
project(umr CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "-O2 -Wall -Wno-array-bounds -fPIC")
endif (CMAKE_BUILD_TYPE STREQUAL "Release")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "-g -O0 -Wall -Wno-array-bounds -DDEBUG -fPIC")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

############
# Packages #
############

find_package(PkgConfig REQUIRED)

find_package(VTK COMPONENTS CommonCore CommonDataModel IOLegacy REQUIRED)

find_package(jsoncpp REQUIRED)
# needed to run on ubuntu (jsoncpp/json, rather than fedora's json)
get_target_property(JSON_INC_PATH jsoncpp_lib INTERFACE_INCLUDE_DIRECTORIES)
include_directories("${JSON_INC_PATH}")

find_package(Boost COMPONENTS program_options unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

#############
# Libraries #
#############

# UMR -- source code library
add_subdirectory(src)
include_directories(src)

# Examples -- examples library
add_subdirectory(examples)

# External libraries (AlgLib)
add_subdirectory(ext)

# Testing
add_subdirectory(tests)

###############
# Executables #
###############

# examples
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/examples)
file(GLOB EX_SOURCES examples/*.cpp)
foreach(filename ${EX_SOURCES})
    get_filename_component(stem ${filename} NAME_WE)
    add_executable(${stem} ${filename})
    target_link_libraries(${stem} PUBLIC UMR AlgLib jsoncpp_lib)
endforeach(filename ${EX_SOURCES})
