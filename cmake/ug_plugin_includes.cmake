# created by Sebastian Reiter

# Plugins are currently build as dynamic libraries.
add_definitions(-DBUILDING_DYNAMIC_LIBRARY)

# the plugin root path
get_filename_component(PLUGIN_ROOT_PATH ${CMAKE_CURRENT_LIST_FILE} PATH)

# include options and dependencies from ug4
set(BUILDING_PLUGIN ON)
include(${PLUGIN_ROOT_PATH}/ug_includes.cmake)

# the plugin shall be written to the plugins folder.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/plugins)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/plugins)

# use uginstall to install the library
# usage: UGInstall("HYPRE" HYPRE_INSTALL_PATH)
# then HYPRE_INSTALL_PATH is the path to hypre.
function(UGInstall name returned_path)
     
   	 # in other cases we use uginstall to download & install HYPRE   	 
	 if(NOT EXISTS "$ENV{UG4_LOCAL_INSTALL_DIR}/${name}/used/installed")
        message(STATUS "Info: ${name} is not provided, but we need ${name}.")
        message(STATUS "Info: Trying to install ${name} with uginstall ${name}.")
        if(NOT defined $ENV{UG4_LOCAL_INSTALL_DIR})
            message(FATAL_ERROR "can't run uginstall. UG4_LOCAL_INSTALL_DIR not available. you need to add 'source YOURUG4DIR/scripts/shell/ugbash' to your .bashrc / .bash_profile")
        endif()
    
        message(STATUS "Install directory is $ENV{UG4_LOCAL_INSTALL_DIR}/${name}.")
        execute_process(COMMAND bash ${UG_ROOT_PATH}/scripts/shell/uginstall ${name})
        if(NOT EXISTS "$ENV{UG4_LOCAL_INSTALL_DIR}/${name}/used/installed")
            message(FATAL_ERROR "${name} could not be installed by uginstall.")
        endif()        
    endif(NOT EXISTS "$ENV{UG4_LOCAL_INSTALL_DIR}/${name}/used/installed")
    
    set(${returned_path} "$ENV{UG4_LOCAL_INSTALL_DIR}/${name}/used/" PARENT_SCOPE)        
endfunction(UGInstall)
