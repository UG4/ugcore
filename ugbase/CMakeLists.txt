cmake_minimum_required(VERSION 2.6)

project(P_UGBASE)

# set global output paths for UGBASE
# They are relative to the source dir for which cmake was executed.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(Sources	ug.cpp)

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase ${Sources})
else(BUILD_ONE_LIB)
	add_library(ugbase ${Sources})
endif(BUILD_ONE_LIB)

########################
# information on the current compilation environment
add_subdirectory(compile_info)

########################
# common
add_subdirectory(common)

########################
# lib_grid
if(buildGrid)
	add_subdirectory(lib_grid)
endif(buildGrid)

########################
# lib_registry
if(buildRegistry)
	add_subdirectory(registry)
endif(buildRegistry)

########################
# bindings
if(buildBindings)
	add_subdirectory(bindings/lua)
	if(buildForVRL)
    	add_subdirectory(bindings/vrl)
	endif(buildForVRL)
endif(buildBindings)

########################
# lib_algebra
if(buildAlgebra)
	add_subdirectory(lib_algebra)	
endif(buildAlgebra)

########################
# lib_disc
if(buildDisc)
	add_subdirectory(lib_disc)
endif(buildDisc)


########################
# pcl
if(PARALLEL)
	add_subdirectory(pcl)
endif(PARALLEL)

########################
# bridge
if(buildBridge)
	add_subdirectory(bridge)
endif(buildBridge)

########################
# ug4 shell
if(buildUGShell)
    add_subdirectory(ug_shell)
endif(buildUGShell)

