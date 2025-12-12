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

#include "demangle.h"

#ifdef UG_POSIX
#include <cxxabi.h>
#endif

#include <sstream>
#include <string>
#include <cstdlib>

using namespace std;

namespace ug {

#ifdef UG_POSIX
string demangle_block(const char *str)
{
	stringstream ss;
	char lastc=0x00;
    string s;
	int status;
	for(char c = *str; c != 0x00; c = *(++str))
	{
		// mangled names start with _ . consider only if last sign was space or tab or newline
		if(c == '_' && (lastc==' ' || lastc == '\n' || lastc == '\t'))
		{
			s = "_";
			// add all characters to the string until space, tab or newline
			for(c = *(++str); c != 0x00; c = *(++str))
			{
				if(c == ' ' || c == '\n' || c == '\t')
					break;
				s += c;
			}
			// some compilers add an additional _ in front. skip it.
			const char *p = s.c_str();
			if(s.length() > 2 && s[1] == '_') p = p+1;
			char *realname = abi::__cxa_demangle(p, 0, 0, &status);
			if(status==0)
				ss << realname; // demangle successfull
			else
				ss << s; // demangling failed, print normal string
			free(realname);
		}
		ss << c;
		lastc =c;
	}
	return ss.str();
}

string demangle(const char *str)
{
	int status;
	char *realname = abi::__cxa_demangle(str, 0, 0, &status);
	if(status==0)
		return realname; // demangle successfull
	else
		return str; // demangling failed, print normal string
}
#else
string demangle_block(const char *str)
{
	return str;
}
string demangle(const char *str)
{
	return str;
}
#endif
}
