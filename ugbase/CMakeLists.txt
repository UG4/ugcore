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
# lib_grid
add_subdirectory(ug_script)

########################
# lib_algebra and lib_discretization

if(NO_ALGEBRA)
	message("-- INFO: Compiling WITHOUT algebra and discretization. (Enable compilation by defining -DNO_ALGEBRA=OFF)")
else(NO_ALGEBRA)
	find_path (HYPRELIB_DIR libHYPRE.a 
			PATHS ENV PATH
			PATH_SUFFIXES src/hypre/lib/ )

	if(HYPRELIB_DIR)
		message("-- HYPRE - Library found at '${HYPRELIB_DIR}'.\n   Compiling WITH Hypre.")

		# add include path
		include_directories (${HYPRELIB_DIR}/../include/)

		# add lib path
		link_directories(${HYPRELIB_LIB_DIR})

		# let all executables link against libHYPRE.a 
		link_libraries (${HYPRELIB_DIR}/libHYPRE.a)

		# add include directories
		include_directories(${P_UGBASE_SOURCE_DIR})

	else(HYPRELIB_DIR)
		message("-- WARNING: HYPRE - Library NOT found. Compiling WITHOUT Hypre.")
	endif(HYPRELIB_DIR)


	if(Boost_FOUND OR HYPRELIB_DIR)
		message("-- INFO: compiling 'algebra' library enabled. (To disable algebra define -DNO_ALGEBRA=ON)")
		add_subdirectory(lib_algebra)	

		message("-- INFO: compiling 'discretization' library enabled")
		add_subdirectory(lib_discretization)
	else(Boost_FOUND OR HYPRELIB_DIR)
		message("-- INFO: Boost-Library and Hypre-Library NOT found. Skipping compilation of 'algebra' and 'discretization'. ")
	endif(Boost_FOUND OR HYPRELIB_DIR)
endif(NO_ALGEBRA)


########################
# pcl

if(MPI_FOUND)
	if(NO_PCL)
		message("-- INFO: not compiling 'pcl' library. (Enable compilation by defining -DNO_PCL=OFF)")
	else(NO_PCL)
		message("-- INFO: compiling 'pcl' library. (Disable compilation by defining -DNO_PCL=ON)")
		add_subdirectory(pcl)
	endif(NO_PCL)

else(MPI_FOUND)
	if(NO_PCL)
		message("-- INFO: not compiling 'pcl' library. (Enable compilation by defining -DNO_PCL=OFF)")
	else(NO_PCL)
		message("-- INFO: not compiling 'pcl' library. MPI is required!")
	endif(NO_PCL)

endif(MPI_FOUND)
