# included from ug_includes.cmake
########################################
# OPENMP
IF(OPENMP)
	# writing to CMAKE_CXX_FLAGS directly can cause problems on some platforms.
	# Please use add_definitions instead. Hope that works in this case, too. sreiter.
  # The linker option '-lgomp' or '-liomp5' can not be added to add_definitions
  # as these are not passed to the link then. But they have to. tklatt.
	#	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")
  IF(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_cxx_flags("-fopenmp")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")
    ADD_DEFINITIONS(-DUG_OPENMP)
    MESSAGE(STATUS "Info: Using OpenMP (experimental)")
  ELSEIF(CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    add_cxx_flags("-fopenmp")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -liomp5")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -liomp5")
    ADD_DEFINITIONS(-DUG_OPENMP)
    MESSAGE(STATUS "Info: Using OpenMP (experimental)")
  ELSEIF(CMAKE_C_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    MESSAGE(WARNING "Clang does not support OpenMP yet.")
  ELSE()
    MESSAGE(WARNING "Don't know compiler type, thus don't know how to enable OpenMP")
  ENDIF()
ENDIF(OPENMP)
