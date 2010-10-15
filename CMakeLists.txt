cmake_minimum_required(VERSION 2.6)



################################################
# ug4
################################################
project(P_UG4)

# set global output paths for UG4
# They are relative to the source dir for which cmake was executed.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# this file is only executed for library builds
# The option triggers a different linking-method in
# "ug_cmake_includes.txt"
set(UG_LIBRARY_BUILD on)

include("ug_cmake_includes.txt")



########################
# ugbase

add_subdirectory(ugbase)

if(NOT UG_SHARED)
    add_subdirectory(apps/ugshell)
endif(NOT UG_SHARED)