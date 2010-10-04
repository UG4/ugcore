# not in use for ug4

cmake_minimum_required(VERSION 2.6)

####
# TetGen Library
####

project(TETGEN)

add_definitions(-DTETLIBRARY)

set_source_files_properties(predicates.cxx
	PROPERTIES COMPILE_FLAGS -O0)

set_source_files_properties(tetgen.cxx
	PROPERTIES COMPILE_FLAGS -Wno-write-strings)

add_library(tet
	tetgen.h
	tetgen.cxx
	predicates.cxx
)
