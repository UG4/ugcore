# Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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
add_definitions(-DUG_TARGET="${TARGET}") # (dummy) preprocessor directive used for displaying build configuration

# The target option enables specific build variables
if("${TARGET}" STREQUAL "ugshell")
	set(buildUGShell ON)
	set(buildForLUA ON)
	set(buildAlgebra ON)
	set(buildPlugins ON)
	set(buildGrid ON)
	set(buildDisc ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)
	set(buildTools ON)
	
elseif("${TARGET}" STREQUAL "libug4")
	set(buildAlgebra ON)
	set(buildPlugins ON)
	set(buildGrid ON)
	set(buildDisc ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)
	set(buildTools ON)

elseif("${TARGET}" STREQUAL "vrl")
	# The vrl works only, if a dynamic library is built.
	if(STATIC_BUILD)
		message(FATAL_ERROR "ug4 for vrl can only be build as a dynamic library. Please set STATIC_BUILD = OFF.")
	endif(STATIC_BUILD)
	
	set(buildAlgebra ON)
	set(buildForVRL ON)
	set(buildForLUA ON)
	set(buildPlugins ON)
	#todo: rename targetLibraryName to ug4_vrl
	#set(targetLibraryName ug4_vrl)
	set(buildGrid ON)
	set(buildDisc ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)
	set(EMBEDDED_PLUGINS ON)
	set(buildTools ON)

elseif("${TARGET}" STREQUAL "vrlgrid")
	if(STATIC_BUILD)
		message(FATAL_ERROR "ug4 for vrl can only be build as a dynamic library. Please set STATIC_BUILD = OFF.")
	endif(STATIC_BUILD)

	set(targetLibraryName vrlgrid)
	set(buildForVRL ON)
	set(buildForLUA ON)
	set(buildGrid ON)
	set(buildRegistry ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildPlugins ON)
	set(BLAS OFF)
	set(LAPACK OFF)
	set(PARALLEL OFF)
	set(DIM 3)
	set(CPU 1)
	set(EMBEDDED_PLUGINS ON)
	set(buildTools ON)

elseif("${TARGET}" STREQUAL "libgrid")
	set(targetLibraryName grid)
	set(buildGrid ON)
	set(buildRegistry ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildPlugins ON)
	set(buildTools ON)
	set(buildForLUA ON)

elseif("${TARGET}" STREQUAL "gridshell")
	set(targetLibraryName ugGrid)
	set(targetExecutableName gridshell)

	set(buildPlugins ON)
	set(buildUGShell ON)
	set(buildForLUA ON)
	set(buildAlgebra OFF)
	set(buildGrid ON)
	set(buildDisc OFF)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)
	set(buildTools ON)

elseif("${TARGET}" STREQUAL "ugplugin")
	set(buildAlgebra ON)
	set(buildGrid ON)
	set(buildDisc ON)
	set(buildForLUA ON)
	set(buildPlugins ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)

elseif("${TARGET}" STREQUAL "amg")
	set(buildPluginSystem OFF)
	set(targetLibraryName ugamg)
	set(buildAlgebra ON)
	set(buildGrid OFF)
	#set(buildPluginSystem OFF)
	set(PARALLEL OFF)
	set(buildTools ON)

elseif("${TARGET}" STREQUAL "luashell")
	set(targetLibraryName ugLua)
	set(targetExecutableName luashell)

	set(buildPlugins OFF)
	set(buildUGShell ON)
	set(buildForLUA ON)
	set(buildAlgebra OFF)
	set(buildGrid OFF)
	set(buildDisc OFF)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)
	set(STATIC_BUILD ON)
	set(EMBEDDED_PLUGINS ON)
	set(DIM 3)
	set(CPU 1)
	set(PARALLEL OFF)
	set(BLAS OFF)
	set(LAPACK OFF)
	set(buildTools ON)

else("${TARGET}" STREQUAL "ugshell")
	message(FATAL_ERROR "Unsupported TARGET: ${TARGET}. Options are: ${targetOptions}")
	
endif("${TARGET}" STREQUAL "ugshell")
