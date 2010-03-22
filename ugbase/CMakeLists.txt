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

#######################
# gcc compiler flags

# Enable these flags for a fast ublas algebra
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -DBOOST_UBLAS_NDEBUG -O3 -Wall")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall")

########################
# common
add_subdirectory(common)

########################
# lib_grid
add_subdirectory(lib_grid)

########################
# lib_algebra and lib_discretization

if(WITH_ALGEBRA)

	find_package(Boost)

	if(Boost_FOUND)
		message("-- Boost - Library found at '${Boost_INCLUDE_DIRS}'. \n   Compiling WITH Boost.")

		# add include path	
		include_directories (${Boost_INCLUDE_DIRS})

	else(Boost_FOUND)
		message("-- WARNING: Boost - Libraray NOT found. Compiling WITHOUT Boost.")
	endif(Boost_FOUND)


	find_path (HYPRELIB_LIB_DIR libHYPRE.a 
			PATHS ENV PATH
			PATH_SUFFIXES src/hypre/lib/ )

	if(HYPRELIB_LIB_DIR)
		message("-- HYPRE - Library found at '${HYPRELIB_LIB_DIR}'.\n   Compiling WITH Hypre.")

		# add include path
		include_directories (${HYPRELIB_LIB_DIR}/../include/)

		# add lib path
		link_directories(${HYPRELIB_LIB_DIR})

		# let all executables link against libHYPRE.a 
		link_libraries (${HYPRELIB_LIB_DIR}/libHYPRE.a)

		# add include directories
		include_directories(${P_UGBASE_SOURCE_DIR})

	else(HYPRELIB_LIB_DIR)
		message("-- WARNING: HYPRE - Libraray NOT found. Compiling WITHOUT Hypre.")
	endif(HYPRELIB_LIB_DIR)

	if(Boost_FOUND OR HYPRE_LIB_DIR)
		message("-- INFO: compiling 'algebra' library enabled")
		add_subdirectory(lib_algebra)	

		message("-- INFO: compiling 'discretization' library enabled")
		add_subdirectory(lib_discretization)
	else(Boost_FOUND OR HYPRE_LIB_DIR)
		message("-- INFO: Boost-Library and Hypre-Library NOT found. Skipping compilation of 'algebra' and 'discretization'. ")
	endif(Boost_FOUND OR HYPRE_LIB_DIR)

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
