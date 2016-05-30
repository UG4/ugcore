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


# Note that you may use -DMPI_DIR=... to set a custom MPI path.

if(PARALLEL)
	add_definitions(-DUG_PARALLEL)

	# search mpi
	if(NOT BUILTIN_MPI)
		find_package(MPI)
		
		if(MPI_FOUND)
			include_directories(${MPI_INCLUDE_PATH})
			# Add mpi libraries:
			# Standard case: add cxx libraries
			if(MPI_CXX_LIBRARIES)
				set(linkLibraries ${linkLibraries} ${MPI_CXX_LIBRARIES})
			# Depreciated case: In order to support cmake versions < 2.8.6, 
			#                   where MPI_CXX_LIBRARIES cannot be used
			else(MPI_CXX_LIBRARIES)				
				if(MPI_LIBRARY)
					set(linkLibraries ${linkLibraries} ${MPI_LIBRARY})
				endif(MPI_LIBRARY)
				if(MPI_EXTRA_LIBRARY)
					set(linkLibraries ${linkLibraries} ${MPI_EXTRA_LIBRARY})
				endif(MPI_EXTRA_LIBRARY)
			endif(MPI_CXX_LIBRARIES)				
		else(MPI_FOUND)
			message(FATAL_ERROR "MPI not found. Please specify a path to your MPI installation "
								"through -DMPI_PATH=... or set PARALLEL to OFF (run cmake -DPARALLEL=OFF .)")
		endif(MPI_FOUND)
	endif(NOT BUILTIN_MPI)

	# MPI is required for parallel builds
endif(PARALLEL)