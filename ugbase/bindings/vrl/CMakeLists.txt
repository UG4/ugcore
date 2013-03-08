cmake_minimum_required(VERSION 2.6)

####
# bindings_vrl Library
####

project(P_BINDINGS_VRL)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(${P_UGBASE_SOURCE_DIR})

include("../../../cmake/ug_includes.cmake")


set(srcBindingsVRL
                    bindings_vrl.cpp
                    type_converter.cpp
                    canvas.cpp
                    messaging.cpp
                    invocation.cpp
                    playground.cpp
                    user_data.cpp
                    basicTest.cpp
)


set(Sources ${srcBindingsVRL})

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/bindings/vrl  ${Sources})
else(BUILD_ONE_LIB)
	add_library(bindings_vrl  ${Sources})
endif(BUILD_ONE_LIB)
