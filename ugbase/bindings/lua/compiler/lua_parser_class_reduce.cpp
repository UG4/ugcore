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

using namespace std;

namespace ug{

nodeType *
LUAParserClass::reduce(nodeType *p)
{
	nodeType *p1, *p2;
	if (!p) return 0;
	switch (p->type)
	{
		case typeCon:
			return p;
			break;
		case typeId:
			return p;
			//cout << id2variable[p->id.i].c_str();
			break;
		case typeOpr:
			switch (p->opr.oper)
			{
				case LUAPARSER_IF:
					p1 = reduce(p->opr.op[0]);
					p2 = reduce(p->opr.op[1]);

					if (p1 != p->opr.op[0])
					{
						freeNode(p->opr.op[0]);
						p->opr.op[0] = p1;
					}
					if (p2 != p->opr.op[1])
					{
						freeNode(p->opr.op[1]);
						p->opr.op[0] = p2;
					}

					break;
				case '=':
					p1 = reduce(p->opr.op[1]);
					if (p1 != p->opr.op[1])
					{
						freeNode(p->opr.op[1]);
						p->opr.op[1] = p1;
					}
					break;

				case 'R':
					p1 = reduce(p->opr.op[0]);
					if (p1 != p->opr.op[0])
					{
						freeNode(p->opr.op[0]);
						p->opr.op[0] = p1;
					}
					break;
				case LUAPARSER_UMINUS:
					p1 = reduce(p->opr.op[0]);
					if (p1->type == typeCon)
					{
						freeNode(p);
						return con(-(p1->con.value));
					}
					if (p1 != p->opr.op[0])
					{
						freeNode(p->opr.op[0]);
						p->opr.op[0] = p1;
					}
					break;

				case LUAPARSER_MATH_COS:
					p1 = reduce(p->opr.op[0]);
					if (p1->type == typeCon)
					{
						freeNode(p);
						return con(cos(p1->con.value));
					}
					if (p1 != p->opr.op[0])
					{
						freeNode(p->opr.op[0]);
						p->opr.op[0] = p1;
					}
					break;

				case LUAPARSER_MATH_SIN:
					p1 = reduce(p->opr.op[0]);
					if (p1->type == typeCon)
					{
						freeNode(p);
						return con(sin(p1->con.value));
					}
					if (p1 != p->opr.op[0])
					{
						freeNode(p->opr.op[0]);
						p->opr.op[0] = p1;
					}
					break;

				case LUAPARSER_MATH_EXP:
					p1 = reduce(p->opr.op[0]);
					if (p1->type == typeCon)
					{
						freeNode(p);
						return con(exp(p1->con.value));
					}
					if (p1 != p->opr.op[0])
					{
						freeNode(p->opr.op[0]);
						p->opr.op[0] = p1;
					}
					break;

				default:
					p1 = reduce(p->opr.op[0]);
					p2 = reduce(p->opr.op[1]);

					if (p->opr.oper != ';' &&
							p1->type == typeCon && p2->type == typeCon)
					{
						double a = p1->con.value;
						double b = p2->con.value;
						freeNode(p);

						switch (p->opr.oper)
						{
							case '+': return con(a + b);
							case '-': return con(a - b);
							case '*': return con(a * b);
							case '/': return con(a / b);
							case LUAPARSER_GE: return con(a > b);
							case LUAPARSER_LE: return con(a < b);
							case LUAPARSER_NE: return con(a != b);
							case LUAPARSER_EQ: return con(a == b);
							case LUAPARSER_AND: return con(a != 0.0 && b != 0.0 ? 1.0: 0.0);
							case LUAPARSER_OR: return con(a != 0.0 || b != 0.0 ? 1.0: 0.0);
						}
					}
					if (p1 != p->opr.op[0])
					{
						freeNode(p->opr.op[0]);
						p->opr.op[0] = p1;
					}
					if (p2 != p->opr.op[1])
					{
						freeNode(p->opr.op[1]);
						p->opr.op[0] = p2;
					}

			}
	}
	return p;
}

}
