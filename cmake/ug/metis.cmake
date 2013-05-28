# included from ug_includes.cmake
########################################
# METIS
if(METIS OR PARMETIS)
	add_definitions(-DUG_METIS)
	include_directories(${UG_ROOT_PATH}/externals/metis-5.0.2/include)
	set(buildMetis ON)
	set(linkLibraries ${linkLibraries} metis)
endif(METIS OR PARMETIS)