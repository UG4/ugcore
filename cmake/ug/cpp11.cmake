# included from ug_includes.cmake
########################################
# C++11
# Note 1: C++11 is still experimental with most compilers, though the most
# important features are already supported.
# Note 2: Though Clang supports C++11, adding the compiler flag globally (as done
# for GCC below) will break the build on trying to compile C/ObjC files with this
# flag. GCC only emits warnings, Clang throws an error.
# TODO: Workaround Note 2. Probably requires larger refactoring of all CMake files.
IF(CXX11)
	IF(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
		# Check for the GCC's C++11 capabilities
		INCLUDE(CheckCXXCompilerFlag)
		CHECK_CXX_COMPILER_FLAG(-std=c++0x HAVE_CXX0X)
		# since GCC4.7 (c++0x will be removed in future versions of GCC)
		CHECK_CXX_COMPILER_FLAG(-std=c++11 HAVE_CXX11)
		# Add appropriate compiler flags
		IF(HAVE_CXX11)
			SET(CXX11_FLAG "-std=c++11")
		ELSEIF(HAVE_CXX0X)
			SET(CXX11_FLAG "-std=c++0x")
		ENDIF()

		IF(CXX11_FLAG)
			ADD_DEFINITIONS(-DUG_CXX11)
			add_cxx_flag(${CXX11_FLAG})
			MESSAGE(STATUS "Info: C++11 enabled. (flag: ${CXX11_FLAG})")
		ELSE()
			SET(CXX11 OFF)
			MESSAGE(STATUS "Info: Compiler does not support C++11 standard.")
		ENDIF()
	ELSE()
		MESSAGE(STATUS "Info: Enabling C++11 is only supported with GCC currently.")
	ENDIF()
ENDIF(CXX11)