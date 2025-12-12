/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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

#include "lua_parser_class.h"
#include "common/assert.h"
#include "lua_compiler_debug.h"
#include "common/profiler/profiler.h"

using namespace std;

namespace ug {

int LUAParserClass::createJITSG(ostream &out, eReturnType r, set<string> &subfunctions)
{
	PROFILE_BEGIN_GROUP(LUAParserClass_createJITSG, "LUA2C");
    int numRet=0;
    switch(r)
    {
        case RT_DIRICHLET: numRet = 1; break;
        case RT_SOURCE: numRet = 1; break;
        case RT_DIFFUSION: numRet = 4; break;
		case RT_VELOCITY: numRet = 2; break;
		case RT_NEUMANN: numRet = 1; break;
        default: UG_ASSERT(0,"");
    }
    if(num_out() != numRet)
    {
    	UG_DLOG(DID_LUACOMPILER, 1, "ERROR in LUA2C for LUA function " << name << ": number of return values must be " << numRet << ", not " << num_out() << "\n");
        return false;
    }
    out << "// INSERTED CODE:";
    out << "\n// JITSG Generated Code from LUA function " << name << "\n";


    for(set<size_t>::iterator it = localFunctions.begin(); it != localFunctions.end(); ++it)
        subfunctions.insert(id2variable[*it]);

    // escape variables
    for(map<size_t, string>::iterator it = id2variable.begin(); it != id2variable.end(); ++it)
		if(is_local((*it).first))
			(*it).second = "_" + (*it).second;
    // TODO: also escape functions

    // rename arguments to x, y, t
    nodeType *a = args;
    id2variable[a->opr.op[0]->id.i] = "x";
    a = a->opr.op[1];
    id2variable[a->opr.op[0]->id.i] = "y";
    id2variable[a->opr.op[1]->id.i] = "t";

    //------ local variables --------
	out << "\t// local variables:\n";
	for(map<size_t, string>::iterator it = id2variable.begin(); it != id2variable.end(); ++it)
		if(is_local((*it).first))
			out << "\tdouble " << (*it).second << ";\n";

    returnType = r;

	out << "\t// code:\n";
	for(size_t i=0; i<nodes.size(); i++)
		createC(nodes[i], out, 1);
	out << "\n// END JITSG Generated Code from LUA function " << name << "\n";
	switch(r)
	{
		case RT_DIRICHLET: out << "dirichletReturn:\n"; break;
        case RT_SOURCE: out << "sourceReturn:\n"; break;
        case RT_DIFFUSION: out << "diffusionReturn:\n"; break;
        case RT_VELOCITY: out << "velocityReturn:\n"; break;
		case RT_NEUMANN: out << "return false;"; break; // necessary?
		default: UG_ASSERT(0,"?");
	}
	out << ";\n";
	return true;
}

}
