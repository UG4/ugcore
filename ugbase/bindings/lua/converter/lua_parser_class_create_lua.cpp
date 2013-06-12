/*
 * lua_parser_class_create_lua.cpp
 *
 *  Created on: 12.06.2013
 *      Author: mrupp
 */
#include "lua_parser_class.h"
#include "common/assert.h"

using namespace std;

namespace ug{

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
				case IF:
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
				case UMINUS:
					out << "-";
					createLUA(p->opr.op[0], out);
					break;

				case MATH_COS:
                case MATH_SIN:
                case MATH_EXP:
                case MATH_ABS:
                case MATH_LOG:
                case MATH_LOG10:
                case MATH_SQRT:
                case MATH_FLOOR:
                case MATH_CEIL:
                    switch (p->opr.oper)
                    {
                        case MATH_COS: out << "math.cos("; break;
                        case MATH_SIN: out << "math.sin("; break;
                        case MATH_EXP: out << "math.exp("; break;
                        case MATH_ABS: out << "math.abs("; break;
                        case MATH_LOG: out << "math.log("; break;
                        case MATH_LOG10: out << "math.log10("; break;
                        case MATH_SQRT: out << "math.sqrt("; break;
                        case MATH_FLOOR: out << "math.floor("; break;
                        case MATH_CEIL: out << "math.ceil("; break;
                    }
					createLUA(p->opr.op[0], out);
					out << ")";
					break;

                case MATH_POW:
                case MATH_MIN:
                case MATH_MAX:
                    switch (p->opr.oper)
                    {
                        case MATH_POW: out << "math.pow("; break;
                        case MATH_MIN: out << "math.min("; break;
                        case MATH_MAX: out << "math.max("; break;
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
						case GE: out << " >= ";
							break;
						case LE: out << " <= ";
							break;
						case NE: out << " ~= ";
							break;
						case EQ: out << " == ";
							break;
						case AND: out << " and ";
							break;
						case OR: out << " or ";
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
