cmake_minimum_required(VERSION 2.6)

####
# uginterface Library
####

project(P_UGINTERFACE)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# This include directory is required by lib_grid.
include_directories(${P_UGBASE_SOURCE_DIR})

set(srcUGInterface	registry.cpp
					interface_base.cpp
					ug_interface.cpp
					standard_interfaces/lib_grid_interface.cpp
					standard_interfaces/lib_discretization_interface.cpp
					ugbridge/class_name_provider.cpp)

# if algebra is enabled, we have to add additional src-files and link with the
# discretization module.
if(NOT NO_ALGEBRA)
	# add the source files
	
else(NOT NO_ALGEBRA)

endif(NOT NO_ALGEBRA)
add_library(uginterface  ${srcUGInterface} ${srcLua})
