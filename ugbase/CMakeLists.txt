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
# lib_algebra and lib_discretization

if(WITH_ALGEBRA)
	message("-- INFO: compiling 'algebra' library enabled")
	add_subdirectory(lib_algebra)	

	message("-- INFO: compiling 'discretization' library enabled")
	add_subdirectory(lib_discretization)
else(WITH_ALGEBRA)
	message("-- INFO: not compiling 'algebra' library (enable option using: 'cmake -DWITH_ALGEBRA=ON')")
	message("-- INFO: not compiling 'discretization' library (enable option using: 'cmake -DWITH_ALGEBRA=ON')")
endif(WITH_ALGEBRA)

########################
# pcl
# pcl will be compiled only if a parallel toolchain has been specified.
if(COMPILE_WITH_PCL)
	# compile pcl
	message("-- INFO: compiling 'pcl' library enabled")
	add_subdirectory(pcl)

else(COMPILE_WITH_PCL)
	# give a warning that pcl will not be compiled.
	message("-- INFO: not compiling 'pcl' library (check the README for more information)")
endif(COMPILE_WITH_PCL)
