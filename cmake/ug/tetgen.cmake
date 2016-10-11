# Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
# Authors: Sebastian Reiter, Martin Rupp
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
if(TETGEN)
	unset(TETGEN_LIBS CACHE)
	find_library(TETGEN_LIBS NAMES tet PATHS ${TETGEN})
	if(TETGEN_LIBS-NOTFOUND)
		message(FATAL_ERROR "ERROR: Couldn't find TETGEN in the specified path.")
	else(TETGEN_LIBS-NOTFOUND)
		add_definitions(-DUG_TETGEN -DTETLIBRARY)
		include_directories(${TETGEN})
		set(linkLibraries ${linkLibraries} ${TETGEN_LIBS})
	endif(TETGEN_LIBS-NOTFOUND)
endif(TETGEN)