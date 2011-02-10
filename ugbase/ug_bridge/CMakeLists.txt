cmake_minimum_required(VERSION 2.6)

####
# ug_bridge Library
####

project(P_UGBRIDGE)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# This include directory is required by lib_grid.
include_directories(${P_UGBASE_SOURCE_DIR})

set(srcUGBridge    	class_name_provider.cpp
					class_helper.cpp
					ug_bridge.cpp
					bridges/lib_grid_bridge.cpp
					bridges/test_bridge.cpp)

# if algebra is enabled, we have to add additional src-files and link with the
# discretization module.
if(NOT NO_ALGEBRA)
	# add source files that depend on the algebra and discretization system here.
	set(srcUGBridge ${srcUGBridge}
					bridges/lib_algebra_bridge.cpp
					bridges/lib_discretization_bridge.cpp
					bridges/user_data_bridge.cpp
					bridges/domain_bridge.cpp)
else(NOT NO_ALGEBRA)

endif(NOT NO_ALGEBRA)

set(Sources ${srcUGBridge})
EXPORTSOURCES(ug_bridge  ${Sources})

if(NOT UG_SHARED)
add_library(ugbridge  ${Sources})
endif(NOT UG_SHARED)
