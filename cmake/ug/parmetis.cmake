# included from ug_includes.cmake
########################################
# PARMETIS
if(PARMETIS)
	add_definitions(-DUG_PARMETIS)
	include_directories(${UG_ROOT_PATH}/externals/parmetis-4.0.2/include)
	set(buildParmetis ON)
	set(linkLibraries ${linkLibraries} parmetis)
endif(PARMETIS)
