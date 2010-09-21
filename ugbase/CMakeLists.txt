################################################
# created by Sebastian Reiter, Andreas Vogel
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
# node_tree
add_subdirectory(node_tree)

########################
# lib_grid
add_subdirectory(lib_grid)

########################
# lib_interface
add_subdirectory(ug_interface)

########################
# lib_script
add_subdirectory(ug_script)

########################
# lib_algebra and lib_discretization

if(NO_ALGEBRA)
	message(STATUS "INFO: Compiling WITHOUT algebra and discretization. (Enable compilation by defining -DNO_ALGEBRA=OFF)")
else(NO_ALGEBRA)
	message(STATUS "INFO: compiling 'algebra' library enabled. (To disable algebra define -DNO_ALGEBRA=ON)")
	add_subdirectory(lib_algebra)	

	message(STATUS "INFO: compiling 'discretization' library enabled")
	add_subdirectory(lib_discretization)
endif(NO_ALGEBRA)


########################
# pcl

if(MPI_FOUND)
	if(NO_PCL)
		message(STATUS "INFO: not compiling 'pcl' library. (Enable compilation by defining -DNO_PCL=OFF)")
	else(NO_PCL)
		message(STATUS "INFO: compiling 'pcl' library. (Disable compilation by defining -DNO_PCL=ON)")
		add_subdirectory(pcl)
	endif(NO_PCL)

else(MPI_FOUND)
	if(NO_PCL)
		message(STATUS "INFO: not compiling 'pcl' library. (Enable compilation by defining -DNO_PCL=OFF)")
	else(NO_PCL)
		message(STATUS "INFO: not compiling 'pcl' library. MPI is required!")
	endif(NO_PCL)

endif(MPI_FOUND)
