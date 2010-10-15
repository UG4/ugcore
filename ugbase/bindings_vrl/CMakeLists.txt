cmake_minimum_required(VERSION 2.6)

####
# bindings_vrl Library
####

project(P_BINDINGS_VRL)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# This include directory is required by lib_grid.
include_directories(${P_UGBASE_SOURCE_DIR})

include("../../ug_cmake_includes.txt")


set(srcBindingsVRL
                    bindings_vrl.cpp
                    type_converter.cpp
                    canvas.cpp
                    messaging.cpp
)


set(Sources ${srcBindingsVRL})
EXPORTSOURCES(bindings_vrl  ${Sources})

if(NOT UG_SHARED)
add_library(bindings_vrl  ${Sources})
endif(NOT UG_SHARED)