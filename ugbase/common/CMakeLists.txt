cmake_minimum_required(VERSION 2.6)

################################################
# common
################################################

project(P_COMMON)

include_directories(${P_UGBASE_SOURCE_DIR})

set(srcCommon	log.cpp
				allocators/small_object_allocator.cpp
                util/string_util.cpp
				util/loader/loader_util.cpp
				util/loader/loader_obj.cpp
				util/ostream_util.cpp
				util/parameter_parsing.cpp
				util/stream_pack.cpp
				math/math_vector_matrix/math_matrix.cpp
				math/math_vector_matrix/math_vector.cpp
				math/misc/tri_box.cpp
				math/misc/tritri.cpp
				math/misc/eigenvalues.cpp
				math/misc/math_util.cpp)
				
set(srcProfiler	profiler/src/ShinyManager.cpp
				profiler/src/ShinyNode.cpp
				profiler/src/ShinyNodePool.cpp
				profiler/src/ShinyOutput.cpp
				profiler/src/ShinyTools.cpp)


EXPORTSOURCES("common" ${srcCommon})
EXPORTSOURCES("common" ${srcProfiler})

if(NOT UG_SHARED)
add_library(common ${srcCommon} ${srcProfiler})
endif(NOT UG_SHARED)