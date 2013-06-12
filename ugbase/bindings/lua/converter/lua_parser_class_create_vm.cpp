/*
 * lua_parser_class_create_vm.cpp
 *
 *  Created on: 12.06.2013
 *      Author: mrupp
 */

#include "lua_parser_class.h"
#include "common/assert.h"

namespace ug{

int LUAParserClass::createVM(nodeType *p, VMAdd &vm)
{
	if (!p) return 0;
	switch (p->type)
	{
		case typeCon:
			vm.push(p->con.value);
			break;
		case typeId:
			vm.push_var(p->id.i);
			break;
		case typeOpr:
			switch (p->opr.oper)
			{
				case IF:
                {
                	// create if directive
                	createVM(p->opr.op[0], vm);

                	//
                	int jmpElsePos = vm.jmp_if_false();

                	std::vector<int> jmpExitPos;
                	createVM(p->opr.op[1], vm);
                	jmpExitPos.push_back( vm.jmp() );

                    nodeType *a = p->opr.op[2];
                    while(a != NULL && a->opr.oper == ELSEIF)
                    {
                    	vm.adjust_jmp_pos(jmpElsePos, vm.get_pos());

                    	// create if directive
                    	createVM(a->opr.op[0], vm);

                    	jmpElsePos = vm.jmp_if_false();

                    	createVM(a->opr.op[1], vm);

                    	jmpExitPos.push_back( vm.jmp() );

                    	a = a->opr.op[2];
                    }
                    if(a != NULL)
                    {
                        UG_ASSERT(a->opr.oper == ELSE, a->opr.oper);
                        vm.adjust_jmp_pos(jmpElsePos, vm.get_pos());

                        createVM(a->opr.op[0], vm);
                    }
                    else
                    	vm.adjust_jmp_pos(jmpElsePos, vm.get_pos());

                    for(size_t i=0; i<jmpExitPos.size(); i++)
                    	vm.adjust_jmp_pos(jmpExitPos[i], vm.get_pos());

                    break;
                }

				case '=':
                    UG_ASSERT(is_local(p->opr.op[0]->id.i), "global variable " << id2variable[p->opr.op[0]->id.i] << " is read-only");

					createVM(p->opr.op[1], vm);

					vm.assign(p->opr.op[0]->id.i);

					break;

                case 'C':
                {
                	UG_ASSERT(0, "not implemented");
                }

				case 'R':
				{
					nodeType *a = p->opr.op[0];
					//i=0;
					while(a->type == typeOpr && a->opr.oper == ',')
					{
						createVM(a->opr.op[0], vm);
						a = a->opr.op[1];
					}
					createVM(a, vm);
					vm.ret();
					break;
				}

				case TK_FOR:
					UG_ASSERT(0, "not implemented");
					break;

				case TK_BREAK:
					UG_ASSERT(0, "not implemented");
					break;

				case UMINUS:
					createVM(p->opr.op[0], vm);
					vm.push(-1.0);
					vm.binary('*');
					break;

                case MATH_PI:
                    vm.push(3.141);
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
                	createVM(p->opr.op[0], vm);
                	vm.unary(p->opr.oper);
					break;

				case '+':
				case '-':
				case '*':
				case '/':
				case '<':
				case '>':
				case GE:
				case LE:
				case NE:
				case EQ:
				case AND:
				case OR:
                case MATH_POW:
                case MATH_MIN:
                case MATH_MAX:
                	createVM(p->opr.op[1], vm);
                	createVM(p->opr.op[0], vm);
                	vm.binary(p->opr.oper);
					break;

                case ';':
					createVM(p->opr.op[0], vm);
					createVM(p->opr.op[1], vm);
					break;
                default:
					UG_ASSERT(0, p->opr.oper  << "not implemented");
					break;
			}
	}
	return 0;
}

}
