/*
 * Copyright (c) 2018:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


/** @file
 * This is a lightweight wrapper for ugshell_main. It simply calls ugshell_main
 * with the provided arguments.
 *
 * If you'd like to do something before ugshell_main is executed, you may
 * write your own main.cpp and compile ug with the
 * \code
 *	-DALTERNATE_MAIN=PATH_TO_YOUR_MAIN/main.cpp
 * \endcode
 * option.
 *
 * Alternatively, you may also create a plugin or enhance your existing one. To this end
 * add an 'alternate_main.cpp' file and a CMakeLists.txt to your plugin which contains the line:
 * \code
 *	set(ALTERNATE_MAIN "${CMAKE_CURRENT_SOURCE_DIR}/alternate_main.cpp" CACHE PATH "Wrapper to ugshell_main." FORCE)
 * \endcode
 *
 * In that case, the alternate_main.cpp is automatically activated as soon as you enable your plugin.
 *
 * If you want to return to the original main, please deactivate the plugin in question and
 * once call:
 * \code
 *	cmake -DUNSET_ALTERNATE_MAIN=ON .
 * \endcode
 */

#include "ug_shell/ugshell_main.h"

int main (int argc, char** argv)
{
	return ugshell_main (argc, argv);
}
