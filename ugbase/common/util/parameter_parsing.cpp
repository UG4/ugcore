/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include "parameter_parsing.h"

#include <cstring>
#include <cstdlib>



namespace ug {

////////////////////////////////////////////////////////////////////////
int GetParamIndex(const char* param, int argc, const char * const * argv) {
	for (int i = 0; i < argc; ++i) {
		if (strcmp(param, argv[i]) == 0) {
			return i;
		}
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
bool FindParam(const char* param, int argc, const char * const * argv) {
	return GetParamIndex(param, argc, argv) != -1;
}


////////////////////////////////////////////////////////////////////////
bool ParamToInt(int& iOut, const char* param, int argc, const char * const * argv) {
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	iOut = atoi(argv[i + 1]);
	return true;
}

////////////////////////////////////////////////////////////////////////
bool ParamToDouble(double &dOut, const char *param, int argc, const char * const * argv)
{
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	dOut = atof(argv[i + 1]);
	return true;
}


////////////////////////////////////////////////////////////////////////
bool ParamToString(const char ** strOut, const char* param, int argc, const char * const * argv) {
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	*strOut = argv[i + 1];
	return true;
}

int ParamToInt(const char* param, int argc, const char * const * argv, int iDefault)
{
	int i;
	if(ParamToInt(i, param, argc, argv) == true)
		return i;
	else return iDefault;
}
double ParamToDouble(const char *param, int argc, const char * const * argv, double dDefault)
{
	double d;
	if(ParamToDouble(d, param, argc, argv) == true)
		return d;
	else return dDefault;
}

}// end of namespace
