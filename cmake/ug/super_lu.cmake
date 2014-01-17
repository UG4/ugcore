
if(UG_CMAKE_SUPER_LU_INCLUDED)
	return()
endif()
set(UG_CMAKE_SUPER_LU_INCLUDED on)

# if super lu is ON
if(SUPERLU)
    set(bla "$ENV{HOME}")
    if(SUPERLU STREQUAL "ON")
        set(SUPERLU_PATH "$ENV{HOME}/local/SuperLU/used/")
    else(SUPERLU STREQUAL "ON")
        set(SUPERLU_PATH SUPERLU)
    endif(SUPERLU STREQUAL "ON")
    message(STATUS "Info: SuperLU is ON, using path ${SUPERLU_PATH}")
    
    add_definitions(-DUG_SUPERLU)
    include_directories("${SUPERLU_PATH}/SRC/")
    set(linkLibraries ${linkLibraries} "${SUPERLU_PATH}/lib/libsuperlu_4.3.a")
endif(SUPERLU)