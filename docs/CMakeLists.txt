# Copyright (c) 2012:  G-CSC, Goethe University Frankfurt
# Author: Torbjörn Klatt
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

MESSAGE(STATUS "Info: Adding make target 'gen_docu' for generating Doxygen documentation.")

# Add the directory with the additional CMake module.
LIST(APPEND CMAKE_MODULE_PATH "../cmake")
# The FindDoxygen CMake module is shipped with CMake since ages; nothing new and shiny.
FIND_PACKAGE(Doxygen QUIET)
# Use the DoxygenTargets CMake module to provide the add_doxygen command
INCLUDE(DoxygenTargets)

# Add the make target while using cmake/DoxygenTargets.doxyfile.in, substitute
# the CMake variables therein and use it as an additional Doxyfile, which overrides
# all Doxygen variables from other Doxyfiles.
ADD_DOXYGEN(
	## the main Doxygen configuration file to be used
	## TODO: maybe make this configurable through a CMAKE argument
	${CMAKE_CURRENT_SOURCE_DIR}/doxy_config_ug4_mathjax.txt
	## where should the documentation be generated to (/html, /latex ...)
	OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ug4
	## the make target name to generate the documentation
	DOC_TARGET gen_docu
	## the directory from where Doxygen should be called (this is docs/)
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	## TODO: add UG's version number to documentation
# 	PROJECT_NUMBER ${UG_VERSION}
	## Uncomment to suppress all Doxygen warnings
	## Leave it commented to write warnings to stderr (WARN_LOGFILE from Doxyfile will be ignored)
# 	NO_WARNINGS
	NO_PDF
	## uncomment QUIET to suppress the verbose Doxygen outout
	## TODO: probably only enable it when DEBUG=ON
	QUIET
)
