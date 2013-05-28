# included from ug_includes.cmake
if(HLIBPRO)
#	find_library(HLIBPRO_LIBS NAMES hpro PATHS ${HLIBPRO})
	find_library(HLIBPRO_LIBS NAMES libhpro.dylib PATH_SUFFIXES ../hlibpro-0.13.6/lib/ ${HLIBPRO}) # hlib built as shared lib via 'scons shared=1 static=0' (26092011ih)
	find_path (HLIBPROLIB_DIR libhpro.dylib
		PATHS ENV PATH
		PATH_SUFFIXES ../hlibpro-0.13.6/lib/ )
#	message(STATUS "INFO: Using HLibPro; content of 'HLIBPRO_LIBS' is '${HLIBPRO_LIBS}', content of 'HLIBPROLIB_DIR' is '${HLIBPROLIB_DIR}'.")
	if(HLIBPRO_LIBS-NOTFOUND)
		message(FATAL_ERROR "ERROR: Couldn't find HLIBPRO in the specified path.")
	else(HLIBPRO_LIBS-NOTFOUND)
#		add_definitions(-DHLIBPROLIB_DIR)
		add_definitions(-DUG_HLIBPRO)
		include_directories(${HLIBPROLIB_DIR}/../include/)
		include_directories(${HLIBPROLIB_DIR}/../src/include/)
		set(linkLibraries ${linkLibraries} ${HLIBPRO_LIBS})
	endif(HLIBPRO_LIBS-NOTFOUND)
endif(HLIBPRO)