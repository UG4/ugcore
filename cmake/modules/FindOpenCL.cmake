# Taken from https://gitorious.org/findopencl/
#
# - Try to find OpenCL
# This module tries to find an OpenCL implementation on your system. It supports
# AMD / ATI, Apple and NVIDIA implementations, but should work, too.
#
# To set manually the paths, define these environment variables:
# OpenCL_INCPATH    - Include path (e.g. OpenCL_INCPATH=/opt/cuda/4.0/cuda/include)
# OpenCL_LIBPATH    - Library path (e.h. OpenCL_LIBPATH=/usr/lib64/nvidia)
#
# Once done this will define
#  OPENCL_FOUND        - system has OpenCL
#  OPENCL_INCLUDE_DIRS  - the OpenCL include directory
#  OPENCL_LIBRARIES    - link these to use OpenCL
#
# WIN32 should work, but is untested

FIND_PACKAGE (PackageHandleStandardArgs)

set (OPENCL_VERSION_STRING "0.1.0")
set (OPENCL_VERSION_MAJOR 0)
set (OPENCL_VERSION_MINOR 1)
set (OPENCL_VERSION_PATCH 0)

if (APPLE)

	FIND_LIBRARY (OPENCL_LIBRARIES OpenCL DOC "OpenCL lib for OSX")
	FIND_PATH (OPENCL_INCLUDE_DIRS OpenCL/cl.h DOC "Include for OpenCL on OSX")
	FIND_PATH (_OPENCL_CPP_INCLUDE_DIRS OpenCL/cl.hpp DOC "Include for OpenCL CPP bindings on OSX")

else (APPLE)

	if (WIN32)

		FIND_PATH (OPENCL_INCLUDE_DIRS CL/cl.h)
		FIND_PATH (_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp)

		# The AMD SDK currently installs both x86 and x86_64 libraries
		# This is only a hack to find out architecture
		if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64" )
			set (OPENCL_LIB_DIR "$ENV{ATISTREAMSDKROOT}/lib/x86_64")
		else (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
			set (OPENCL_LIB_DIR "$ENV{ATISTREAMSDKROOT}/lib/x86")
		endif ()
		FIND_LIBRARY (OPENCL_LIBRARIES OpenCL.lib PATHS ${OPENCL_LIB_DIR} ENV OpenCL_LIBPATH)

		GET_FILENAME_COMPONENT (_OPENCL_INC_CAND ${OPENCL_LIB_DIR}/../../include ABSOLUTE)

		# On Win32 search relative to the library
		FIND_PATH (OPENCL_INCLUDE_DIRS CL/cl.h PATHS "${_OPENCL_INC_CAND}" ENV OpenCL_INCPATH)
		FIND_PATH (_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp PATHS "${_OPENCL_INC_CAND}" ENV OpenCL_INCPATH)

	else ()

		# Unix style platforms
		FIND_LIBRARY (OPENCL_LIBRARIES OpenCL
			PATHS ENV LD_LIBRARY_PATH ENV OpenCL_LIBPATH
		)

		GET_FILENAME_COMPONENT (OPENCL_LIB_DIR ${OPENCL_LIBRARIES} PATH)
		GET_FILENAME_COMPONENT (_OPENCL_INC_CAND ${OPENCL_LIB_DIR}/../../include ABSOLUTE)

		# The AMD SDK currently does not place its headers
		# in /usr/include, therefore also search relative
		# to the library
		FIND_PATH (OPENCL_INCLUDE_DIRS CL/cl.h PATHS ${_OPENCL_INC_CAND} "/usr/local/cuda/include" "/opt/AMDAPP/include" ENV OpenCL_INCPATH)
		FIND_PATH (_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp PATHS ${_OPENCL_INC_CAND} "/usr/local/cuda/include" "/opt/AMDAPP/include" ENV OpenCL_INCPATH)

	endif ()

endif ()

FIND_PACKAGE_HANDLE_STANDARD_ARGS (OpenCL DEFAULT_MSG OPENCL_LIBRARIES OPENCL_INCLUDE_DIRS)

if (_OPENCL_CPP_INCLUDE_DIRS)
	set ( OPENCL_HAS_CPP_BINDINGS TRUE )
	LIST ( APPEND OPENCL_INCLUDE_DIRS ${_OPENCL_CPP_INCLUDE_DIRS} )
	# This is often the same, so clean up
	LIST ( REMOVE_DUPLICATES OPENCL_INCLUDE_DIRS )
endif ()

MARK_AS_ADVANCED (
  OPENCL_INCLUDE_DIRS
)
