# Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
# Author: Stephan Grein
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

################################################################################
# Declare a method that allows all sub-cmake-files to add their own linker flags
# to the main project (P_UG4) defines
#######################
# Export flags to global variable.
# PURPOSE: use this function to add local project linker flags to global
#          ug4LibLinkerFlags property which is used to build libug4 project
# @param flags   list of linker flags to be exported
function(ExportLinkerFlags flags)
    # iterate over all arguments and insert given prefix
    foreach(flag ${ARGV})
            # retrieve the global property ugLinkerFlags and store its values
            # in a temp variable
            # NOTE: properties must be assigned to variables before being used
            get_property(temp GLOBAL PROPERTY ugLinkerFlags)
            # append tmp to the global ugDependencies property using the correct prefix
            set_property(GLOBAL PROPERTY ugLinkerFlags "${temp} ${flag}")
    endforeach(flag)
endfunction(ExportLinkerFlags)
