# included from ug_includes.cmake
########################################
# HYPRE
if(HYPRE)
	find_library(HYPRE_LIBS HYPRE PATHS ${HYPRE})

	if(HYPRE_LIBS-NOTFOUND)
		message(FATAL_ERROR "ERROR: Couldn't find HYPRE in the specified path.")
	else(HYPRE_LIBS-NOTFOUND)
		add_definitions(-DUG_HYPRE -DHYPRELIB_DIR) # TODO: Is '-DHYPRELIB_DIR' used?
		include_directories(${HYPRE}/../include/)
		set(linkLibraries ${linkLibraries} ${HYPRE_LIBS})
	endif(HYPRE_LIBS-NOTFOUND)
endif(HYPRE)
