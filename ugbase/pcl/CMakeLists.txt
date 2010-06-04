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
			pcl_methods.cpp)

add_library(pcl ${srcPcl})
