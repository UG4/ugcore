# included from ug_includes.cmake
if(PARALLEL)
	if(BUILTIN_MPI)
		add_definitions(-DUG_PARALLEL)
	else(BUILTIN_MPI)
		# search mpi
		find_package(MPI)
		# MPI is required for parallel builds
		if(MPI_FOUND)
			add_definitions(-DUG_PARALLEL)
			include_directories(${MPI_INCLUDE_PATH})
			# Add mpi libraries:
			# Standard case: add cxx libraries
			if(MPI_CXX_LIBRARIES)
				set(linkLibraries ${linkLibraries} ${MPI_CXX_LIBRARIES})
			# Depreciated case: In order to support cmake versions < 2.8.6, 
			#                   where MPI_CXX_LIBRARIES cannot be used
			else(MPI_CXX_LIBRARIES)				
				set(linkLibraries ${linkLibraries} ${MPI_LIBRARY})
				if(MPI_EXTRA_LIBRARY)
				set(linkLibraries ${linkLibraries} ${MPI_EXTRA_LIBRARY})
				endif(MPI_EXTRA_LIBRARY)
			endif(MPI_CXX_LIBRARIES)				
		else(MPI_FOUND)
			message(FATAL_ERROR "MPI not found. Please set PARALLEL to OFF (run cmake -DPARALLEL=OFF ...)")
		endif(MPI_FOUND)
	endif(BUILTIN_MPI)
endif(PARALLEL)