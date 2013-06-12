/*
 * lua_parser_class_create_jitsg.cpp
 *
 *  Created on: 12.06.2013
 *      Author: mrupp
 */
#include "lua_parser_class.h"
#include "common/assert.h"

using namespace std;

namespace ug{

int LUAParserClass::createJITSG(ostream &out, eReturnType r, set<string> &subfunctions)
{
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
        UG_LOG("ERROR in LUA2C for LUA function " << name << ": number of return values must be " << numRet << ", not " << num_out() << "\n");
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
