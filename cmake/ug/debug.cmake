# included from ug_includes.cmake
message(STATUS "Debug = ${DEBUG}, ${CMAKE_BUILD_TYPE}")
if(DEBUG)
	add_definitions(-DUG_DEBUG)
	# if no build type set, add default debug flags
	if(NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "None" OR
		"${CMAKE_BUILD_TYPE}" STREQUAL "")
		# use cmake standard flags for debug builds (-g)
		add_cxx_flag(${CMAKE_CXX_FLAGS_DEBUG})
		# if user specified additional debug flags add them
		if(DEBUG_FORMAT)
			add_cxx_flag(${DEBUG_FORMAT})
			message(STATUS "Info: Debug Information is ${DEBUG_FORMAT}") 
		endif(DEBUG_FORMAT)
	elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
		message(WARNING "Build type set to Release, but DEBUG build wanted.")
	elseif("${CMAKE_BUILD_TYPE}" STREQUAL "MinSizeRel")
		message(WARNING "Build type set to MinSizeRel, but DEBUG build wanted.")
	endif()
	
	# This code would enable strict bounds checking for STL objects like in vector::operator[]. 
	# however, GLIBCXX_DEBUG and strstream don't work together on mac (bug: http://bit.ly/cH78bC). 
	# when this bug is fixed, one could set those flags (or similar) depending on the compiler.
	
	# I disabled the following definitions, since they cause several additional problems.
	# First they also cause problems on MinGW builds.
	# Secondly they have to be defined also for all executables, which link against
	# ug. Since they are well hidden, this may cause quite some headache.
	# One should think about a flag DEBUG_STL or something like this, which explicitly
	# enables those flags. (sreiter)
	#if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND NOT APPLE)
	#	add_definitions(-D_GLIBCXX_DEBUG=1 -D_GLIBCXX_DEBUG_PEDANTIC=1)
	#ENDIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND NOT APPLE)
else(DEBUG)
	# add release definitions
	add_definitions(-DBOOST_UBLAS_NDEBUG)
	if(NOT CMAKE_BUILD_TYPE)
		# if DEBUG=OFF also add standard cmake release flags (-O3 -DNDEBUG)
		add_cxx_flag(${CMAKE_CXX_FLAGS_RELEASE})
	elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
		message(WARNING "Build type set to Debug, but DEBUG=OFF. Leads to strange cflags!")
	endif()

	# compiler specific release flags
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		add_cxx_flag("-funroll-loops -ftree-vectorize")
	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")
		add_cxx_flag("-hipa5 -hunroll2")
	endif()
endif(DEBUG)

if(DEBUG_LOGS)
	add_definitions(-DUG_ENABLE_DEBUG_LOGS)
endif(DEBUG_LOGS)