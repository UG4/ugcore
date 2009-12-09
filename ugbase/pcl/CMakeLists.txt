################################################
# created by Sebastian Reiter
# s.b.reiter@googlemail.com
#
# Parallel Communication Layer
################################################

cmake_minimum_required(VERSION 2.6)

project(P_PCL)

include_directories(${P_UGBASE_SOURCE_DIR})

add_library(pcl pcl_methods.cpp)
