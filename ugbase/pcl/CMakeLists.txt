################################################
# created by Sebastian Reiter
# s.b.reiter@googlemail.com
#
# Parallel Communication Layer
################################################

cmake_minimum_required(VERSION 2.6)

project(P_PCL)

include_directories(${P_UGBASE_SOURCE_DIR})
include_directories("externals/metis")

set(srcPcl	pcl_base.cpp
			pcl_methods.cpp
			pcl_process_communicator.cpp
			pcl_util.cpp
			parallel_archive.cpp)

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/pcl ${srcPcl})
else(BUILD_ONE_LIB)
	add_library(pcl ${srcPcl})
endif(BUILD_ONE_LIB)
