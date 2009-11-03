################################################
# created by Sebastian Reiter
# s.b.reiter@googlemail.com
#
# ugbase
################################################

cmake_minimum_required(VERSION 2.6)

project(P_UGBASE)

# set global output paths for UGBASE
# They are relative to the source dir for which cmake was executed.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

########################
# common
add_subdirectory(common)

########################
# lib_grid
add_subdirectory(lib_grid)

########################
# pcl
# pcl will be compiled only if a parallel toolchain has been specified.
if(COMPILE_WITH_PCL)
	# compile pcl
	add_subdirectory(pcl)

else(COMPILE_WITH_PCL)
	# give a warning that pcl will not be compiled.
	message("")
	message("WARNING: parallel applications and libraries are not compiled ...")
	message("     ... consider to specify a ParallelToolchain file.\n")
endif(COMPILE_WITH_PCL)
