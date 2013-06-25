/*
 * lua_parser_class_create_vm.cpp
 *
 *  Created on: 12.06.2013
 *      Author: mrupp
 */

#include "lua_parser_class.h"
#include "common/assert.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/info_commands.h"


using namespace std;
namespace ug{

int LUAParserClass::createVM(nodeType *p, VMAdd &vm,  std::map<std::string, SmartPtr<VMAdd> > &subVM)
{
	if (!p) return 0;
	switch (p->type)
	{
		case typeCon:
			vm.push(p->con.value);
			break;
		case typeId:
			if(is_global(p->id.i))
			{
				int i=p->id.i;
				if(id2variable[i].compare("true")==0)
					vm.push(1.0);
				else if(id2variable[i].compare("false")==0)
					vm.push(0.0);
				else
				{
					lua_State* L = ug::script::GetDefaultLuaState();
					vm.push(ug::bridge::LuaGetNumber(L, id2variable[i].c_str(), 0) );
				}
			}
			else
				vm.push_var(p->id.i);
			break;
		case typeOpr:
			switch (p->opr.oper)
			{
				case IF:
                {
                	// create if directive
                	createVM(p->opr.op[0], vm, subVM);

                	//
                	int jmpElsePos = vm.jmp_if_false();

                	std::vector<int> jmpExitPos;
                	createVM(p->opr.op[1], vm, subVM);
                	jmpExitPos.push_back( vm.jmp() );

                    nodeType *a = p->opr.op[2];
                    while(a != NULL && a->opr.oper == ELSEIF)
                    {
                    	vm.adjust_jmp_pos(jmpElsePos, vm.get_pos());

                    	// create if directive
                    	createVM(a->opr.op[0], vm, subVM);

                    	jmpElsePos = vm.jmp_if_false();

                    	createVM(a->opr.op[1], vm, subVM);

                    	jmpExitPos.push_back( vm.jmp() );

                    	a = a->opr.op[2];
                    }
                    if(a != NULL)
                    {
                        UG_ASSERT(a->opr.oper == ELSE, a->opr.oper);
                        vm.adjust_jmp_pos(jmpElsePos, vm.get_pos());

                        createVM(a->opr.op[0], vm, subVM);
                    }
                    else
                    	vm.adjust_jmp_pos(jmpElsePos, vm.get_pos());

                    for(size_t i=0; i<jmpExitPos.size(); i++)
                    	vm.adjust_jmp_pos(jmpExitPos[i], vm.get_pos());

                    break;
                }

				case '=':
                    UG_ASSERT(is_local(p->opr.op[0]->id.i), "global variable " << id2variable[p->opr.op[0]->id.i] << " is read-only");

					createVM(p->opr.op[1], vm, subVM);

					vm.assign(p->opr.op[0]->id.i);

					break;

                case 'C':
                {
                	nodeType *a = p->opr.op[1];
					while(a->type == typeOpr && a->opr.oper == ',')
					{
						createVM(a->opr.op[0], vm, subVM);
						a = a->opr.op[1];
					}
					createVM(a, vm, subVM);

					std::string subName = id2variable[p->opr.op[0]->id.i];
					SmartPtr<VMAdd> sub = subVM[subName];
					UG_ASSERT(sub.valid(), "subroutine " << subName << " not found.");
                	vm.call(sub);
                	break;
                }

				case 'R':
				{
					nodeType *a = p->opr.op[0];
					//i=0;
					while(a->type == typeOpr && a->opr.oper == ',')
					{
						createVM(a->opr.op[0], vm, subVM);
						a = a->opr.op[1];
					}
					createVM(a, vm, subVM);
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
					createVM(p->opr.op[0], vm, subVM);
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
                	createVM(p->opr.op[0], vm, subVM);
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
                	createVM(p->opr.op[1], vm, subVM);
                	createVM(p->opr.op[0], vm, subVM);
                	vm.binary(p->opr.oper);
					break;

                case ';':
					createVM(p->opr.op[0], vm, subVM);
					createVM(p->opr.op[1], vm, subVM);
					break;
                default:
					UG_ASSERT(0, p->opr.oper  << "not implemented");
					break;
			}
	}
	return 0;
}

int LUAParserClass::createVM(VMAdd &vm)
{
	vm.set_name(name);
	std::map<std::string, SmartPtr<LUAParserClass> > subfunctions;
	std::map<std::string, SmartPtr<VMAdd> > subVM;
	add_subfunctions(subfunctions);
	UG_LOG("Subfunctions:\n");
	for(std::map<std::string, SmartPtr<LUAParserClass> >::iterator it = subfunctions.begin(); it != subfunctions.end(); ++it)
	{
		std::string subName = (*it).first;
		UG_LOG(subName << "\n");
		SmartPtr<VMAdd> sub = new VMAdd;
		sub->set_name(subName);
		subVM[subName] = sub;
	}

	for(std::map<std::string, SmartPtr<LUAParserClass> >::iterator it = subfunctions.begin(); it != subfunctions.end(); ++it)
	{
		UG_LOG((*it).first << "\n");
		std::string name = (*it).first;
		((*it).second)->createVMHeader(*subVM[name]);
	}

	createVMHeader(vm);


	for(std::map<std::string, SmartPtr<LUAParserClass> >::iterator it = subfunctions.begin(); it != subfunctions.end(); ++it)
	{
		UG_LOG((*it).first << "\n");
		std::string name = (*it).first;
		((*it).second)->createVMSub(*subVM[name], subVM);
	}


	createVMSub(vm, subVM);
	return true;
}

int LUAParserClass::createVMHeader(VMAdd &vm)
{
	nodeType *a = args;
	int i=0;
	while(a->type == typeOpr)
	{
		UG_LOG(id2variable[a->opr.op[0]->id.i] << " - " << a->opr.op[0]->id.i << "\n");
		a = a->opr.op[1];
		i++;
	}
	UG_LOG(id2variable[a->id.i] << " - " << a->id.i << "\n");

	vm.set_in_out(i+1,numOut);
	vm.set_nr_of_variables(variables.size());
	UG_LOG("function " << name << " in: " << i+1 << " out " << numOut << "\n");

	return true;
}

int LUAParserClass::createVMSub(VMAdd &vm, std::map<std::string, SmartPtr<VMAdd> > &subVM)
{
	UG_LOG("CODE:\n");
	for(size_t i=0; i<nodes.size(); i++)
		createVM(nodes[i], vm, subVM);
	return true;
}

int LUAParserClass::add_subfunctions(std::map<std::string, SmartPtr<LUAParserClass> > &subfunctions)
{
    for(set<size_t>::iterator it = localFunctions.begin(); it != localFunctions.end(); ++it)
        if(add_subfunction(id2variable[*it], subfunctions) == false) return false;
    return true;
}

int LUAParserClass::add_subfunction(std::string name, std::map<std::string, SmartPtr<LUAParserClass> > &subfunctions)
{
    //UG_LOG("adding " << name << "\n");
    if(subfunctions.find(name) != subfunctions.end()) return true;

    SmartPtr<LUAParserClass> parser = new LUAParserClass;

    if(parser->parse_luaFunction(name.c_str()) == false)
        return false;

    if(parser->num_out() != 1)
    {
        UG_LOG("ERROR in LUA2C for LUA function " << name << ":  subfunction must have exactly one return value (not " << parser->num_out() << ")\n");
        return false;
    }

    parser->returnType = RT_SUBFUNCTION;
    subfunctions[name] = parser;

    parser->add_subfunctions(subfunctions);
	return true;
}



}
