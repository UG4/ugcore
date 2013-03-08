# created by Sebastian Reiter

# Plugins are currently build as dynamic libraries.
add_definitions(-DBUILDING_DYNAMIC_LIBRARY)

# the plugin root path
get_filename_component(PLUGIN_ROOT_PATH ${CMAKE_CURRENT_LIST_FILE} PATH)

# include options and dependencies from ug4
set(BUILDING_PLUGIN ON)
include(${PLUGIN_ROOT_PATH}/ug_includes.cmake)

# Plugins only work in a shared environment
if(NOT USE_NEW_CMAKE_INCLUDES)
	if(STATIC)
		message(FATAL_ERROR "Plugins require build option -DUG_SHARED=ON. Aborting.")
	endif(STATIC)
	
	if(NOT UG_PLUGINS)
		message(FATAL_ERROR "Plugins require build option -DUG_PLUGINS=ON. Aborting.")
	endif(NOT UG_PLUGINS)
endif(NOT USE_NEW_CMAKE_INCLUDES)

# the plugin shall be written to the plugins folder.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/plugins)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/plugins)
