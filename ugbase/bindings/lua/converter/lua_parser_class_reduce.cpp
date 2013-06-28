/*
 * lua_parser_class_reduce.cpp
 *
 *  Created on: 12.06.2013
 *      Author: mrupp
 */

#include "lua_parser_class.h"
#include "common/assert.h"

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
