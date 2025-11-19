/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef DEMANGLE_H_
#define DEMANGLE_H_
#include <string>

namespace ug{
/**
 * demangles C++ function names like _ZZ12ug_backtracev = ug_backtrace().
 * also demangles them when a lot of them appear "in between". make sure they are
 * seperated by ' ', '\n' or '\t' and start with _, like they do in backtrace_symbols.
 * Works only in POSIX, otherwise returns str
 * @sa demangle
 * @param str mangled strings, e.g. _ZZ12ug_backtracev
 * @return the demangled string e.g. ug_backtrace()
 */
std::string demangle_block(const char *str);

/**
 * demangles C++ function and class names
 * Works only in POSIX, otherwise returns str
 * @param str mangled string, containing stuff e.g. 3barI5emptyLi17EE
 * @return the demangled string e.g. bar<empty, 17>
 * */
std::string demangle(const char *str);
}
#endif
