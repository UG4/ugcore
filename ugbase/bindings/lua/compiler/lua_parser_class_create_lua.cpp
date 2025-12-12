/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

using namespace std;

namespace ug {

int LUAParserClass::createLUA(nodeType *p, ostream &out)
{
	if (!p) return 0;
    nodeType *a;
    //int i;
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
					out << "if ";
					createLUA(p->opr.op[0], out);
					out << " then\n";
					createLUA(p->opr.op[1], out);
					out << "end\n";

					break;
				case '=':
                    UG_ASSERT(is_local(p->opr.op[0]->id.i), "global variable " << id2variable[p->opr.op[0]->id.i] << " is read-only");
					out << id2variable[p->opr.op[0]->id.i] << " = ";
					createLUA(p->opr.op[1], out);
					out << "\n";
					break;
				case 'R':
                    out << "return ";
                    a = p->opr.op[0];
					//i=0;
					while(a->type == typeOpr && a->opr.oper == ',')
					{
						createLUA(a->opr.op[0], out);
						out << ", ";
						a = a->opr.op[1];
					}
					createLUA(a, out);
					out << "\n";
					break;

                case 'C':
                    out << id2variable[p->opr.op[0]->id.i].c_str() << "(";
                    a = p->opr.op[1];
					//i=0;
					while(a->type == typeOpr && a->opr.oper == ',')
					{
						createLUA(a->opr.op[0], out);
						out << ", ";
						a = a->opr.op[1];
					}
					createLUA(a, out);
					out << ")\n";
					break;
				case LUAPARSER_UMINUS:
					out << "-";
					createLUA(p->opr.op[0], out);
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
                        case LUAPARSER_MATH_COS: out << "math.cos("; break;
                        case LUAPARSER_MATH_SIN: out << "math.sin("; break;
                        case LUAPARSER_MATH_EXP: out << "math.exp("; break;
                        case LUAPARSER_MATH_ABS: out << "math.abs("; break;
                        case LUAPARSER_MATH_LOG: out << "math.log("; break;
                        case LUAPARSER_MATH_LOG10: out << "math.log10("; break;
                        case LUAPARSER_MATH_SQRT: out << "math.sqrt("; break;
                        case LUAPARSER_MATH_FLOOR: out << "math.floor("; break;
                        case LUAPARSER_MATH_CEIL: out << "math.ceil("; break;
                    }
					createLUA(p->opr.op[0], out);
					out << ")";
					break;

                case LUAPARSER_MATH_POW:
                case LUAPARSER_MATH_MIN:
                case LUAPARSER_MATH_MAX:
                    switch (p->opr.oper)
                    {
                        case LUAPARSER_MATH_POW: out << "math.pow("; break;
                        case LUAPARSER_MATH_MIN: out << "math.min("; break;
                        case LUAPARSER_MATH_MAX: out << "math.max("; break;
                    }

					createLUA(p->opr.op[0], out);
					out << ", ";
                    createLUA(p->opr.op[1], out);
                    out << ")";
					break;

				case ';':
					createLUA(p->opr.op[0], out);
					createLUA(p->opr.op[1], out);
					break;

				default:
					out << "(";
					createLUA(p->opr.op[0], out);
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
						case LUAPARSER_NE: out << " ~= ";
							break;
						case LUAPARSER_EQ: out << " == ";
							break;
						case LUAPARSER_AND: out << " and ";
							break;
						case LUAPARSER_OR: out << " or ";
							break;
					}
					out << "(";
					createLUA(p->opr.op[1], out);
					out << ")";
			}
	}
	return 0;
}


int LUAParserClass::createLUA(ostream &out)
{
    out << "function " << name << "(";
    nodeType *a = args;
	while(a->type == typeOpr)
	{
		out << id2variable[a->opr.op[0]->id.i] << ", ";
		a = a->opr.op[1];
	}
	out << id2variable[a->id.i] << ")\n";

	for(size_t i=0; i<nodes.size(); i++)
		createLUA(nodes[i], out);

    out << "end\n";
    return true;
}

}
