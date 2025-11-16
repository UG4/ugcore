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
#include "common/util/string_util.h"
#include "lua_compiler_debug.h"

using namespace std;
namespace ug{

int LUAParserClass::createC(nodeType *p, ostream &out, int indent)
{
	if (!p) return 0;
	switch (p->type)
	{
		case typeCon:
			out << p->con.value;
			break;
		case typeId:
			getVar(p->id.i, out);
			break;
		case typeOpr:
			switch (p->opr.oper)
			{
				case LUAPARSER_IF:
                {
					out << repeat('\t', indent);			out << "if(";
					createC(p->opr.op[0], out, 0);
					out << ")\n";
					out << repeat('\t', indent); 			out << "{\n";
					createC(p->opr.op[1], out, indent+1);
					out << repeat('\t', indent); 			out << "}\n";

                    nodeType *a = p->opr.op[2];
                    while(a != nullptr && a->opr.oper == LUAPARSER_ELSEIF)
                    {
                        out << repeat('\t', indent); 			out << "else if(";
                        createC(a->opr.op[0], out, 0);
                        out << ")\n";
                        out << repeat('\t', indent); 			out << "{\n";
                        createC(a->opr.op[1], out, indent+1);
                        out << repeat('\t', indent); 			out << "}\n";
                        a = a->opr.op[2];
                    }
                    if(a != nullptr)
                    {
                        UG_ASSERT(a->opr.oper == LUAPARSER_ELSE, a->opr.oper);
                        out << repeat('\t', indent); 			out << "else\n";
                        out << repeat('\t', indent); 			out << "{\n";
                        createC(a->opr.op[0], out, indent+1);
                        out << repeat('\t', indent); 			out << "}\n";
                    }

                    break;
                }

				case '=':
                    UG_ASSERT(is_local(p->opr.op[0]->id.i), "global variable " << id2variable[p->opr.op[0]->id.i] << " is read-only");
					out << repeat('\t', indent);
					out << id2variable[p->opr.op[0]->id.i] << " = ";
					createC(p->opr.op[1], out, 0);
					out << ";\n";
					break;

                case 'C':
                {
                    out << "LUA2C_Subfunction_" << id2variable[p->opr.op[0]->id.i] << "(";
                    nodeType *a = p->opr.op[1];
					while(a->type == typeOpr && a->opr.oper == ',')
					{
						createC(a->opr.op[0], out, indent);
						out << ", ";
						a = a->opr.op[1];
					}
					createC(a, out, indent);
					out << ")";
                    break;
                }

				case 'R':
				{
					if(returnType == RT_SUBFUNCTION || returnType == RT_NEUMANN)
                    {
                        nodeType *a = p->opr.op[0];
                        if(a->type == typeOpr && a->opr.oper == ',')
                        {
                            UG_ASSERT(0, "subfunctions may not return more then one value");
                            return false;
                        }
                        out << repeat('\t', indent);
                        out << "return ";
                        createC(a, out, indent);
                        out << ";\n";
                    }
                    else if(returnType == RT_DIFFUSION)
                    {
                        const char *rt[] = {"A11", "A12", "A21", "A22"};
                        createRT(p->opr.op[0], out, rt, 4, indent);

                        out << "\n";
                        out << repeat('\t', indent);
                        out << "goto diffusionReturn;\n";
                    }
                    else if(returnType == RT_VELOCITY)
                    {
                        const char *rt[] = {"vx", "vy"};
                        createRT(p->opr.op[0], out, rt, 2, indent);

                        out << "\n";
                        out << repeat('\t', indent);
                        out << "goto velocityReturn;\n";
                    }
                    else if(returnType == RT_DIRICHLET || returnType == RT_SOURCE)
                    {
                        const char *rt[] = {"f"};
                        createRT(p->opr.op[0], out, rt, 1, indent);

                        out << "\n";
                        out << repeat('\t', indent);
						if(returnType == RT_DIRICHLET)
							out << "goto dirichletReturn;\n";
						else
							out << "goto sourceReturn;\n";
                    }
                    else
                    {
                        nodeType *a = p->opr.op[0];
                        int i=0;
                        while(a->type == typeOpr && a->opr.oper == ',')
                        {
                            out << repeat('\t', indent);
                            out << "LUA2C_ret[" << i++ << "] = ";
                            createC(a->opr.op[0], out, indent);
                            out << ";\n";
                            a = a->opr.op[1];
                        }
                        out << repeat('\t', indent);
                        out << "LUA2C_ret[" << i++ << "] = ";
                        createC(a, out, indent);
                        out << ";\n";
                        out << repeat('\t', indent);
                        out << "return 1;\n";
                    }
                    break;
				}

				case LUAPARSER_FOR:
					out << repeat('\t', indent);
					out << "for(";
					createC(p->opr.op[0], out, 0);
					out << " = ";
					createC(p->opr.op[1], out, 0);
					out << "; ";
					createC(p->opr.op[0], out, 0);
					out << " <= ";
					createC(p->opr.op[2], out, 0);
					out << "; ";
					createC(p->opr.op[0], out, 0);
					out << " += ";
					createC(p->opr.op[3], out, 0);
					out << ")\n";
					out << repeat('\t', indent);
					out << "{\n";
						createC(p->opr.op[4], out, indent+1);
					out << repeat('\t', indent);
					out << "}\n";
					break;

				case LUAPARSER_BREAK:
					out << repeat('\t', indent);
					out << "break;\n";
					break;

				case LUAPARSER_UMINUS:
					out << "-(";
					createC(p->opr.op[0], out, 0);
					out << ")";
					break;

                case LUAPARSER_MATH_PI:
                    out << " LUAPARSER_MATH_PI ";
                    break;

				case LUAPARSER_MATH_COS:
                case LUAPARSER_MATH_SIN:
                case LUAPARSER_MATH_EXP:
                case LUAPARSER_MATH_ABS:
                case LUAPARSER_MATH_LOG:
                case LUAPARSER_MATH_LOG10:
                case LUAPARSER_MATH_SQRT:
                case LUAPARSER_MATH_FLOOR:
                case LUAPARSER_MATH_CEIL:
                    switch (p->opr.oper)
                    {
                        case LUAPARSER_MATH_COS: out << "cos("; break;
                        case LUAPARSER_MATH_SIN: out << "sin("; break;
                        case LUAPARSER_MATH_EXP: out << "exp("; break;
                        case LUAPARSER_MATH_ABS: out << "fabs("; break;
                        case LUAPARSER_MATH_LOG: out << "log("; break;
                        case LUAPARSER_MATH_LOG10: out << "log10("; break;
                        case LUAPARSER_MATH_SQRT: out << "sqrt("; break;
                        case LUAPARSER_MATH_FLOOR: out << "floor("; break;
                        case LUAPARSER_MATH_CEIL: out << "ceil("; break;
                    }
					createC(p->opr.op[0], out, 0);
					out << ")";
					break;

                case LUAPARSER_MATH_POW:
                case LUAPARSER_MATH_MIN:
                case LUAPARSER_MATH_MAX:
                    switch (p->opr.oper)
                    {
                        case LUAPARSER_MATH_POW: out << "pow("; break;
                        case LUAPARSER_MATH_MIN: out << "min("; break;
                        case LUAPARSER_MATH_MAX: out << "max("; break;
                    }

					createC(p->opr.op[0], out, 0);
					out << ", ";
                    createC(p->opr.op[1], out, 0);
                    out << ")";
					break;

				case ';':
					createC(p->opr.op[0], out, indent);
					createC(p->opr.op[1], out, indent);
					break;

				default:
					out << "(";
					createC(p->opr.op[0], out, 0);
					out << ")";

					switch (p->opr.oper)
					{
						case '+': out << '+';
							break;
						case '-': out << '-';
							break;
						case '*': out << '*';
							break;
						case '/': out << '/';
							break;
						case '<': out << '<';
							break;
						case '>': out << '>';
							break;
						case LUAPARSER_GE: out << " >= ";
							break;
						case LUAPARSER_LE: out << " <= ";
							break;
						case LUAPARSER_NE: out << " != ";
							break;
						case LUAPARSER_EQ: out << " == ";
							break;
						case LUAPARSER_AND: out << " && ";
							break;
						case LUAPARSER_OR: out << " || ";
							break;
					}
					out << "(";
					createC(p->opr.op[1], out, 0);
					out << ")";
			}
	}
	return 0;
}

int LUAParserClass::createC(ostream &out)
{
     // local functions

    set<string> knownFunctions;
    stringstream declarations;
    stringstream definitions;
    int ret;
    ret = add_subfunctions(knownFunctions, declarations, definitions);
    if(ret != LUAParserOK)
    {
    	UG_DLOG(DID_LUACOMPILER, 2, "add_subfunctions failed.\n");
        return ret;
    }

    out << "#define LUAPARSER_MATH_PI 3.1415926535897932384626433832795028841971693\n";
    out << "// inline function declarations\n";
    out << declarations.str() << "\n";

    out << "// inline function definitions\n";
    out << definitions.str() << "\n";

    // the function
	out << "int " << name << "(";
	out << "double *LUA2C_ret, "; //[" << numOut << "], ";
	out << "const double *LUA2C_in)\n"; //[" << numIn << "])\n";
	nodeType *a = args;

	out << "{\n";

    int i=0;
	while(a->type == typeOpr)
	{
		out << "\tdouble " << id2variable[a->opr.op[0]->id.i] << " = LUA2C_in[" << i++ << "];\n";
		a = a->opr.op[1];
	}
	out << "\tdouble " << id2variable[a->id.i] << " = LUA2C_in[" << i++ << "];\n";

	//------ local variables --------
	out << "\t// local variables:\n";
	for(map<string, size_t>::iterator it = variables.begin(); it != variables.end(); ++it)
		if(is_local((*it).second))
			out << "\tdouble " << (*it).first << ";\n";


	out << "\t// code:\n";
	for(size_t i=0; i<nodes.size(); i++)
		createC(nodes[i], out, 1);
	out << "}\n";
	return LUAParserOK;
}


int LUAParserClass::createC_inline(ostream &out)
{
    out << "inline ";
    declare(out);
    out << "\n{\n";

	//------ local variables --------
	out << "\t// local variables:\n";
	for(map<string, size_t>::iterator it = variables.begin(); it != variables.end(); ++it)
		if(is_local((*it).second))
			out << "\tdouble " << (*it).first << ";\n";

	out << "\t// code:\n";
	for(size_t i=0; i<nodes.size(); i++)
		createC(nodes[i], out, 1);
	out << "}\n";
	return true;
}

int LUAParserClass::addfunctionC(string name, set<string> &knownFunctions, stringstream &declarations, stringstream &definitions)
{
    //UG_LOG("adding " << name << "\n");
    if(knownFunctions.find(name) != knownFunctions.end()) return true;
    knownFunctions.insert(name);

    LUAParserClass parser;
    int ret = parser.parse_luaFunction(name.c_str());
    if(ret == LUAParserError)
    {
    	UG_DLOG(DID_LUACOMPILER, 1, "ERROR in LUA2C for LUA function " << name << ":  could not be parsed.\n");
        return ret;
    }
    else if(ret != LUAParserOK)
    {
    	return ret;
    }

    if(parser.num_out() != 1)
    {
    	UG_DLOG(DID_LUACOMPILER, 1, "ERROR in LUA2C for LUA function " << name << ":  subfunctions must have exactly one return value (not " << parser.num_out() << ")\n");
        return LUAParserError;
    }

    parser.returnType = RT_SUBFUNCTION;

    declarations << "inline ";
    parser.declare(declarations);
    declarations << ";\n";

    parser.createC_inline(definitions);

    ret = parser.add_subfunctions(knownFunctions, declarations, definitions);
    if(ret != LUAParserOK)
    {
    	UG_DLOG(DID_LUACOMPILER, 1, "ERROR in LUA2C for LUA function " << name << ":  one of its subfunctions could not be parsed.\n");
    	return ret;
    }

    return LUAParserOK;
}


int LUAParserClass::declare(ostream &out)
{
    out << "double LUA2C_Subfunction_" << name << "(";
	nodeType *a = args;
	while(a->type == typeOpr)
	{
		out << "double " << id2variable[a->opr.op[0]->id.i] << ", ";
		a = a->opr.op[1];
	}
	out << "double " << id2variable[a->id.i] << ")";
    return true;
}

}
