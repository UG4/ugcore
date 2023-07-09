# Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
# Authors: Arne Naegel, Tobias Trautmann
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


# This file is used for cpack configuration.
set(CPACK_PACKAGE_NAME "UG4")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/UG4")
set(CPACK_PACKAGE_VENDOR "TheUG4Group")
set(CPACK_PACKAGE_CONTACT "ug4@uni-frankfurt.de")
set(CPACK_PACKAGE_CHECKSUM SHA256)

# make DESTDIR=/home/xyz install
set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/ug4")

#include(cmake/cpack/zip-config.cmake) # Create ZIP (generic)
#include(cmake/cpack/nuget-config.cmake)
#include(cmake/cpack/osx-config.cmake)
#include(cmake/cpack/rpm-config.cmake)
#include(cmake/cpack/deb-config.cmake)

# We can select a generator using cpack -G
# set(CPACK_GENERATOR "DEB") # Override? # ZIP, DEB
include(CPack)

# These are the 'component' categories used (extend, if neccessary!)
cpack_add_component(applications DISPLAY_NAME "Executables")
cpack_add_component(libraries DISPLAY_NAME "Libraries")
cpack_add_component(sources DISPLAY_NAME "Sources")
