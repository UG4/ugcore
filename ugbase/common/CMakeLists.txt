cmake_minimum_required(VERSION 2.6)

################################################
# common
################################################

project(P_COMMON)

include_directories(${P_UGBASE_SOURCE_DIR})

set(srcCommon	log.cpp
				util/loader/loader_util.cpp
				util/loader/loader_obj.cpp
				util/stream_pack.cpp
				math/math_vector_matrix/math_matrix.cpp
				math/math_vector_matrix/math_vector.cpp)
				
set(srcProfiler	profiler/src/ShinyManager.cpp
				profiler/src/ShinyNode.cpp
				profiler/src/ShinyNodePool.cpp
				profiler/src/ShinyOutput.cpp
				profiler/src/ShinyTools.cpp)

add_library(common ${srcCommon} ${srcProfiler})
