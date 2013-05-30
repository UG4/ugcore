# included from ug_includes.cmake
if(CUDA)
	SET(MY_HOSTNAME $ENV{HOSTNAME})
	if(MY_HOSTNAME STREQUAL "tales")
        SET(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda-5.0)
    endif()    
    find_package(CUDA)
	if(CUDA_FOUND)
		SET(CUDA_NVCC_FLAGS "-arch=sm_13")
		MESSAGE(STATUS "Info: CUDA ${CUDA_VERSION}. Toolkit root dir: ${CUDA_TOOLKIT_ROOT_DIR}, nvcc: ${CUDA_NVCC_EXECUTABLE}")

		file(WRITE ${CMAKE_BINARY_DIR}/cuda_jit.h
			"#define CUDA_TOOLKIT_ROOT_DIR \"${CUDA_TOOLKIT_ROOT_DIR}\"\n"
			"#define CUDA_NVCC_EXECUTABLE \"${CUDA_NVCC_EXECUTABLE}\"\n"
			"#define CUDA_VERSION \"${CUDA_VERSION}\"\n"
			"#define CUDA_NVCC_FLAGS \"${CUDA_NVCC_FLAGS}\"\n")

		include_directories(${CUDA_TOOLKIT_INCLUDE})
		include_directories(${CUDA_TOOLKIT_ROOT_DIR}/samples/common/inc)
        
		set(linkLibraries ${linkLibraries} ${CUDA_CUDART_LIBRARY} ${CUDA_cublas_LIBRARY})
		## cusparse is not found on some systems
		#set(linkLibraries ${linkLibraries} "${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcusparse.so")
		set(linkLibraries ${linkLibraries} ${CUDA_cusparse_LIBRARY})   
		add_definitions(-DCUDA_AVAILABLE)
	else(CUDA_FOUND)
		message(FATAL_ERROR "Error: Couldn't find CUDA")
	endif()
endif()