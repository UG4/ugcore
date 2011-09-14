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


include("ug_cmake_includes.txt")



########################
# ugbase

add_subdirectory(ugbase)

# compile core plugins
if(UG_PLUGINS)
	add_subdirectory(plugins/core)
endif(UG_PLUGINS)
