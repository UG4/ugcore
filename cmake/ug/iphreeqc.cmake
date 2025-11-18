# Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
# Author: Ingo Heppner
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# Begin: 23172015ih
# included from ug_includes.cmake
########################################

# IPhreeqc - Modules Based on the Geochemical Model PHREEQC for Use in Scripting and Programming Languages

if (PHREEQC)
	message (STATUS "    ###### Begin of iphreeqc.cmake: ################################################")
#	message(STATUS "    INFO in iphreeqc.cmake: 'PHREEQC'             = '${PHREEQC}' (TMP)")
	
	set (PHREEQC_ROOT_PATH ${UG_ROOT_PATH}externals/iphreeqc-3.3.2-10335)
#	message (STATUS "    INFO in iphreeqc.cmake: 'UG_ROOT_PATH'        = '${UG_ROOT_PATH}' (TMP)")
#	message (STATUS "    INFO in iphreeqc.cmake: 'PHREEQC_ROOT_PATH'   = '${PHREEQC_ROOT_PATH}' - TMP)")
	
	if (DEBUG)
		set (PHREEQCLIB_BASENAME IPhreeqcd) # base name of IPhreeqc lib if built with CMAKE_BUILD_TYPE = 'Debug'
	else ()
		set (PHREEQCLIB_BASENAME IPhreeqc)  # base name of IPhreeqc lib if built with CMAKE_BUILD_TYPE = ''
	endif ()
#	message(STATUS "    INFO in iphreeqc.cmake: 'PHREEQCLIB_BASENAME' = '${PHREEQCLIB_BASENAME}' - TMP)")

	# Force 'find_library()' to find the static version of the IPhreeqC library (cf.
	# http://stackoverflow.com/questions/16344302/cmake-ignores-static-library-link-request;
	# for the Windows related settings cf.
	# http://stackoverflow.com/questions/3762057/cmake-how-to-produce-binaries-as-static-as-possible.
	# NOTE: After switching to or from '-DSTATIC_BUILD=ON' you might have to delete your CMakeCache.txt):
#	message (STATUS "    INFO in iphreeqc.cmake: CMAKE_FIND_LIBRARY_SUFFIXES before resetting: " ${CMAKE_FIND_LIBRARY_SUFFIXES} )
	if (STATIC_BUILD)
		if (WIN32)
			set (CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
		else()
			set (CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
		endif ()
	endif ()

	message (STATUS "    INFO in iphreeqc.cmake: CMAKE_FIND_LIBRARY_SUFFIXES after  resetting: " ${CMAKE_FIND_LIBRARY_SUFFIXES} )
	
	# Find library:
	find_library (PHREEQC_LIBS NAMES ${PHREEQCLIB_BASENAME} PATHS ${PHREEQC_ROOT_PATH}/build/)
	message (STATUS "    INFO in iphreeqc.cmake: After  executing 'find_library()':        'PHREEQC_LIBS'        = '${PHREEQC_LIBS}' (TMP)")
	if (PHREEQC_LIBS-NOTFOUND)
		message (FATAL_ERROR "    ERROR in iphreeqc.cmake: Couldn't find PHREEQC in the specified path.")

	else()
		add_definitions (-DUG_PHREEQC)
		get_property (inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#		message(STATUS "    INFO in iphreeqc.cmake: Before executing 'include_directories()': 'INCLUDE_DIRECTORIES' = '${inc_dirs}'; TMP)")
		include_directories (${PHREEQC_ROOT_PATH}/src
			            ${PHREEQC_ROOT_PATH}/src/phreeqcpp/common
				    ${PHREEQC_ROOT_PATH}/src/phreeqcpp/PhreeqcKeywords)
		get_property (inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#		message (STATUS "    INFO in iphreeqc.cmake: After  executing 'include_directories()': 'INCLUDE_DIRECTORIES' = '${inc_dirs}'; TMP)")
#		message (STATUS "    INFO in iphreeqc.cmake: Before setting   'linkLibraries': 'linkLibraries' = '${linkLibraries}' (TMP)")
		set (linkLibraries ${linkLibraries} ${PHREEQC_LIBS})
#		message (STATUS "    INFO in iphreeqc.cmake: After  setting   'linkLibraries': 'linkLibraries' = '${linkLibraries}' (TMP)")
	endif ()
	message (STATUS "    ###### End   of iphreeqc.cmake: ################################################")
endif ()
