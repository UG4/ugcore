# Copyright (c) 2011-2014:  G-CSC, Goethe University Frankfurt
# Author: Sebastian Reiter
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
        if(NOT DEFINED ENV{UG4_LOCAL_INSTALL_DIR})
            message(FATAL_ERROR "can't run uginstall. UG4_LOCAL_INSTALL_DIR not available. you need to add 'source YOURUG4DIR/scripts/shell/ugbash' to your .bashrc / .bash_profile."
            " if you already added it, try source .bashrc")
        endif()
    
        message(STATUS "Install directory is $ENV{UG4_LOCAL_INSTALL_DIR}/${name}.")
        execute_process(COMMAND bash ${UG_ROOT_PATH}/ugcore/scripts/shell/uginstall ${name})
        if(NOT EXISTS "$ENV{UG4_LOCAL_INSTALL_DIR}/${name}/used/installed")
            message(FATAL_ERROR "${name} could not be installed by uginstall.")
        endif()        
    endif(NOT EXISTS "$ENV{UG4_LOCAL_INSTALL_DIR}/${name}/used/installed")
    
    set(${returned_path} "$ENV{UG4_LOCAL_INSTALL_DIR}/${name}/used/" PARENT_SCOPE)        
endfunction(UGInstall)
