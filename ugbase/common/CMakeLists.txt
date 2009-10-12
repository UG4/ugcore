cmake_minimum_required(VERSION 2.6)

################################################
# common
################################################

project(P_COMMON)

set(srcCommon	common.cpp
				util/loader/loader_util.cpp
				util/loader/loader_obj.cpp
				math/math_vector_matrix/math_matrix.cpp
				math/math_vector_matrix/math_vector.cpp)

add_library(common ${srcCommon})
