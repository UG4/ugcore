cmake_minimum_required(VERSION 2.6)

####
# registry Library
####

project(P_REGISTRY)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include("../../cmake/ug_includes.cmake")

set(srcRegistry    	registry.cpp
                    class.cpp
					class_name_provider.cpp
					global_function.cpp
					class_helper.cpp
					registry_util.cpp)

set(Sources ${srcRegistry})

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/registry  ${Sources})
else(BUILD_ONE_LIB)
	add_library(registry  ${Sources})
endif(BUILD_ONE_LIB)
