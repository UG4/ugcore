# Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
# Author: Markus Breit
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

##########################################################
# Allows all sub-cmake-files to add their own flags to   #
# specified files when building with EMBEDDED_PLUGINS=ON #
##########################################################

# USAGE
# variant 1:
# exportSingleFileCompileFlags(file "def1;def2;def3;...")
# 
# variant 2:
# set(defList def1 def2 def3 ...)
# exportSingleFileCompileFlags(file "${defList}")
#
# example: exportSingleFileCompileFlags(cpp11_feature.cpp "-std=c++11;-DSOME_PREPROCESSOR_FLAG")
#
# file must be a valid file, prepended by its path, typically ${CMAKE_CURRENT_SOURCE_DIR};
#
# @param file the file definitions are to be exported for
# @param definitions a list of definitions to be exported
function(exportSingleFileCompileFlags file flags)
	list(LENGTH flags len)
	if (${len} GREATER 0)
		math(EXPR len "${len} - 1")
		foreach (idx RANGE ${len})
			list(GET flags ${idx} flag)
			set_property(GLOBAL APPEND PROPERTY ugSingleFileCompileFlagsFiles ${file})
			set_property(GLOBAL APPEND PROPERTY ugSingleFileCompileFlagsDefs "${flag}")
		endforeach()
	endif()
endfunction(exportSingleFileCompileFlags)
