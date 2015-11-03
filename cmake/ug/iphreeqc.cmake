# Begin: 23172015ih
# included from ug_includes.cmake
########################################

# IPhreeqc - Modules Based on the Geochemical Model PHREEQC for Use in Scripting and Programming Languages

if(PHREEQC)
	message(STATUS "    ###### Begin of iphreeqc.cmake: ################################################")
#	message(STATUS "    INFO in iphreeqc.cmake: 'PHREEQC'             = '${PHREEQC}' (TMP)")
	
	set(PHREEQC_ROOT_PATH ${UG_ROOT_PATH}externals/iphreeqc-3.3.2-10335)
#	message(STATUS "    INFO in iphreeqc.cmake: 'UG_ROOT_PATH'        = '${UG_ROOT_PATH}' (TMP)")
#	message(STATUS "    INFO in iphreeqc.cmake: 'PHREEQC_ROOT_PATH'   = '${PHREEQC_ROOT_PATH}' - TMP)")
	
	if(DEBUG)
		set(PHREEQCLIB_BASENAME IPhreeqcd) # base name of IPhreeqc lib if built with CMAKE_BUILD_TYPE = 'Debug'
	else(DEBUG)
		set(PHREEQCLIB_BASENAME IPhreeqc)  # base name of IPhreeqc lib if built with CMAKE_BUILD_TYPE = ''
	endif(DEBUG)
#	message(STATUS "    INFO in iphreeqc.cmake: 'PHREEQCLIB_BASENAME' = '${PHREEQCLIB_BASENAME}' - TMP)")

	# Force 'find_library()' to find the static version of the IPhreeqC library (cf.
	# http://stackoverflow.com/questions/16344302/cmake-ignores-static-library-link-request;
	# for the Windows related settings cf.
	# http://stackoverflow.com/questions/3762057/cmake-how-to-produce-binaries-as-static-as-possible.
	# NOTE: After switching to or from '-DSTATIC_BUILD=ON' you might have to delete your CMakeCache.txt):
#	MESSAGE(STATUS "    INFO in iphreeqc.cmake: CMAKE_FIND_LIBRARY_SUFFIXES before resetting: " ${CMAKE_FIND_LIBRARY_SUFFIXES} )
	if (STATIC_BUILD)
		IF(WIN32)
			SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
		ELSE(WIN32)
			SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
		ENDIF(WIN32)
	endif(STATIC_BUILD)
	MESSAGE(STATUS "    INFO in iphreeqc.cmake: CMAKE_FIND_LIBRARY_SUFFIXES after  resetting: " ${CMAKE_FIND_LIBRARY_SUFFIXES} )
	
	# Find library:
	find_library(PHREEQC_LIBS NAMES ${PHREEQCLIB_BASENAME} PATHS ${PHREEQC_ROOT_PATH}/build/)
	message(STATUS "    INFO in iphreeqc.cmake: After  executing 'find_library()':        'PHREEQC_LIBS'        = '${PHREEQC_LIBS}' (TMP)")
	if(PHREEQC_LIBS-NOTFOUND)
		message(FATAL_ERROR "    ERROR in iphreeqc.cmake: Couldn't find PHREEQC in the specified path.")
	else(PHREEQC_LIBS-NOTFOUND)
		add_definitions(-DUG_PHREEQC)
		get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#		message(STATUS "    INFO in iphreeqc.cmake: Before executing 'include_directories()': 'INCLUDE_DIRECTORIES' = '${inc_dirs}'; TMP)")
		include_directories(${PHREEQC_ROOT_PATH}/src
			            ${PHREEQC_ROOT_PATH}/src/phreeqcpp/common
				    ${PHREEQC_ROOT_PATH}/src/phreeqcpp/PhreeqcKeywords)
		get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#		message(STATUS "    INFO in iphreeqc.cmake: After  executing 'include_directories()': 'INCLUDE_DIRECTORIES' = '${inc_dirs}'; TMP)")
#		message(STATUS "    INFO in iphreeqc.cmake: Before setting   'linkLibraries': 'linkLibraries' = '${linkLibraries}' (TMP)")
		set(linkLibraries ${linkLibraries} ${PHREEQC_LIBS})
#		message(STATUS "    INFO in iphreeqc.cmake: After  setting   'linkLibraries': 'linkLibraries' = '${linkLibraries}' (TMP)")
	endif(PHREEQC_LIBS-NOTFOUND)
	message(STATUS "    ###### End   of iphreeqc.cmake: ################################################")
endif(PHREEQC)
