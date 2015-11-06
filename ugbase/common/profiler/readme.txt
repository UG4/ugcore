# Copyright (c) 2009:  G-CSC, Goethe University Frankfurt
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

The profiler source code has been taken from the shiny-profiler project by Aidin Abedi.
See http://shinyprofiler.sourceforge.net/ and have a look at the "Shiny  C++ Profiler.htm" in this folder.
Some files have been slightly altered to ensure compatibility with ug and to remove compilation errors.
altered files are:
	ShinyMacros.h
	ShinyPrereqs.h
	ShinyNodePool.cpp
	ShinyTools.h

ShinyTools.h defines a hash-function. The original implementation which simply casted the pointer to and uint32_t did not compile on 64-bit machines. I added a small workaround. Sadly this will most likely introduce errors in the hashing on 64 bit platforms.
If you are using a 64 bit platform and experience strange profiler-behaviour you should definitively think about improving hashing in shiny.

The profiler can be disabled by adding
#define SHINY_PROFILER FALSE
before you include profiler.h

Sebastian	(s.b.reiter@googlemail.com)
