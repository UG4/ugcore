# created by Sebastian Reiter
# s.b.reiter@googlemail.com
# y09 m26 d06

cmake_minimum_required(VERSION 2.6)

################################################################
# grid-unit-tests
################################################################

project(P_GRID_UNIT_TESTS)

# set global output paths for the sample_app
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# This include directory is required by lib_grid.
include_directories(${P_GRID_UNIT_TESTS_SOURCE_DIR}/../ugbase)

# The path in which ug's libs are located - if they have already been built.
link_directories(${P_GRID_UNIT_TESTS_SOURCE_DIR}/../lib)

# create the executable
add_executable(gut src/main.cpp)

# the order in which the libraries are specified is cruical!
# we're linking with libgrid and libcommon (libgrid depends on libcommon).
target_link_libraries(gut grid common)
