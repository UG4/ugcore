cmake_minimum_required(VERSION 2.6)

####
# ug_bridge Library
####

project(P_UGBRIDGE)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include("../../ug_cmake_includes.txt")

set(srcUGBridge    	ug_bridge.cpp
					lib_grid_bridge.cpp
					test_bridge.cpp
					profiler_bridge.cpp
					util_bridge.cpp
					misc_bridge.cpp
					domain_bridge.cpp)

set(srcUGBridge     ${srcUGBridge}
					pcl_bridge.cpp)

# if algebra is enabled, we have to add additional src-files and link with the
# discretization module.
if(NOT NO_ALGEBRA)
	# add source files that depend on the algebra and discretization system here.
	set(srcUGBridge ${srcUGBridge}
					lib_algebra_bridges/lib_algebra_bridge.cpp
					lib_algebra_bridges/amg_bridge.cpp
					
					lib_disc_bridges/lib_disc_common.cpp
					lib_disc_bridges/lib_disc_algebra.cpp
					lib_disc_bridges/lib_disc_domain.cpp
					lib_disc_bridges/lib_disc_disc.cpp
					
					user_data_bridge.cpp)
else(NOT NO_ALGEBRA)

endif(NOT NO_ALGEBRA)

set(Sources ${srcUGBridge})
EXPORTSOURCES(ug_bridge  ${Sources})

if(NOT UG_SHARED)
add_library(ugbridge  ${Sources})
endif(NOT UG_SHARED)
