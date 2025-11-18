# Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
# Author: Martin Rupp
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

# included from ug_includes.cmake
if (HLIBPRO)
#	find_library(HLIBPRO_LIBS NAMES hpro PATHS ${HLIBPRO})
	find_library (HLIBPRO_LIBS NAMES libhpro.dylib PATH_SUFFIXES ../hlibpro-0.13.6/lib/ ${HLIBPRO}) # hlib built as shared lib via 'scons shared=1 static=0' (26092011ih)
	find_path (HLIBPROLIB_DIR libhpro.dylib
		PATHS ENV PATH
		PATH_SUFFIXES ../hlibpro-0.13.6/lib/ )
#	message (STATUS "INFO: Using HLibPro; content of 'HLIBPRO_LIBS' is '${HLIBPRO_LIBS}', content of 'HLIBPROLIB_DIR' is '${HLIBPROLIB_DIR}'.")

	if (HLIBPRO_LIBS-NOTFOUND)
		message (FATAL_ERROR "ERROR: Couldn't find HLIBPRO in the specified path.")

	else ()
#		add_definitions(-DHLIBPROLIB_DIR)
		add_definitions (-DUG_HLIBPRO)
		include_directories (${HLIBPROLIB_DIR}/../include/)
		include_directories (${HLIBPROLIB_DIR}/../src/include/)
		set (linkLibraries ${linkLibraries} ${HLIBPRO_LIBS})
	endif ()
endif ()