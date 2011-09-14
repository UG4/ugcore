cmake_minimum_required(VERSION 2.6)

################################################
# common
################################################

project(P_COMMON)

include("../../ug_cmake_includes.txt")

set(srcCommon	log.cpp
				serialization.cpp
				allocators/small_object_allocator.cpp
				util/binary_buffer.cpp
                util/string_util.cpp
				util/loader/loader_util.cpp
				util/loader/loader_obj.cpp
				util/ostream_util.cpp
				util/parameter_parsing.cpp
				math/math_vector_matrix/math_matrix.cpp
				math/math_vector_matrix/math_vector.cpp
				math/misc/tri_box.cpp
				math/misc/tritri.cpp
				math/misc/eigenvalues.cpp
				math/misc/math_util.cpp)

set(sources ${srcCommon})


if(UG_PROFILER)
	set(srcProfiler	profiler/profiler.cpp
					profiler/src/ShinyManager.cpp
					profiler/src/ShinyNode.cpp
					profiler/src/ShinyNodePool.cpp
					profiler/src/ShinyOutput.cpp
					profiler/src/ShinyTools.cpp)
	
	set(sources ${sources} ${srcProfiler})
endif(UG_PROFILER)


################################################################################
# Platform dependend code
################################################################################
if(UNIX)
	set(sources ${sources} os_dependent/file_util_posix.cpp)
	
	if(MAC)
		set(sources ${sources} os_dependent/os_info_osx.cpp)
	else(MAC)
		set(sources ${sources} os_dependent/os_info_linux.cpp)
	endif(MAC)

	if(UG_PLUGINS)
		set(sources ${sources} os_dependent/plugin_util_unix.cpp)
	endif(UG_PLUGINS)

elseif(WIN32)
	set(sources ${sources} os_dependent/os_info_win.cpp)
	
	if(MINGW)
		set(sources ${sources} os_dependent/file_util_posix.cpp)
		
	else(MINGW)
	endif(MINGW)
	
	if(UG_PLUGINS)
		message(FATAL_ERROR "Plugins not yet supported for MSVC. Aborting.")
	endif(UG_PLUGINS)
	
	
else(UNIX)
	message(STATUS "YOUR OS IS NOT FULLY SUPPORTED. Methods in common/os_dependent won't be available!")
	if(UG_PLUGINS)
		message(FATAL_ERROR "Plugins not yet supported for this platform. Consider deactivating plugins through cmake -DUG_PLUGINS=OFF")
	endif(UG_PLUGINS)
endif(UNIX)


################################################################################
#	add sources to library
################################################################################
EXPORTSOURCES("common" ${sources})

if(NOT UG_SHARED)
	add_library(common ${sources})
endif(NOT UG_SHARED)