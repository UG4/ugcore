/* 
 * File:   lua_parser_class.h
 * Author: mrupp
 *
 * Created on 20. November 2012, 10:16
 */

#ifndef LUAParserClass_H
#define	LUAParserClass_H

#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <string>
#include <math.h>

#include "parser_node.h"
#include "parser.hpp"
#include <iostream>
#include "common/util/smart_pointer.h"

#include "vm.h"

void yyerror(const char *s);
namespace ug{
class LUAParserClass;
}
void yaccparse(const char*command, ug::LUAParserClass *p);

namespace ug{


class LUAParserClass
{
	std::map<std::string, size_t> variables;
	std::map<size_t, std::string> id2variable;
	std::set<size_t> localVariables;
    std::set<size_t> localFunctions;

	std::vector<nodeType *> nodes;
	std::string name;
	int numOut;
	nodeType *args;

public:    
    enum eReturnType
    {
        RT_SUBFUNCTION, RT_DIFFUSION, RT_VELOCITY, RT_CALLBACK, RT_DIRICHLET, 
        RT_SOURCE, RT_NEUMANN
    };
    eReturnType returnType;
public:
	LUAParserClass()
	{
        numOut = -1;
        returnType = RT_CALLBACK;
		args = NULL;
	}
	
	void set_name(const char *s)
	{
		name = s;
	}
	
	void set_name(int id)
	{
		set_name(get_name_for_id(id));
	}
	
	
	
	void set_arguments(nodeType *p)
	{
		args = p;
	}
	
	std::stringstream err;
	
	void
	add(nodeType *p)
	{
		nodes.push_back(p);
	}
	
	void parse(const char *command)
	{
		clear();
		yaccparse(command, this);
	}

	void
	clear()
	{
		variables.clear();
		id2variable.clear();
		localVariables.clear();
		err.clear();
		nodes.clear();
	}

	const char *
	get_name_for_id(size_t id)
	{
		return id2variable[id].c_str();
	}

	int
	get_id_for_name(const char*name);
	
	bool
	set_local(size_t id)
	{
		localVariables.insert(id);
		return true;
	}
	
    bool
    is_global(size_t id)
    {
        return !is_local(id) && !is_arg(id);
    }

	bool
	is_local(size_t id)
	{
		return localVariables.find(id) != localVariables.end();
	}
	
	bool
	is_arg(size_t id)
	{
		nodeType *a = args;	
		while(a->type == typeOpr)
		{
			if(a->opr.op[0]->id.i == (int)id) return true;
			a = a->opr.op[1];
		}
		if(a->id.i == (int)id) return true;
		return false;
	}


	void
	assert_local(size_t id)
	{
		if (!is_local(id))
			err << "Error: " << get_name_for_id(id) << " not local! Try adding \"local " << get_name_for_id(id) << "\"\n";
	}

	bool
	has_errors()
	{
		return err.str().length() != 0;
	}
	
	int num_in()
	{
		nodeType *a = args;
		int i=0;
		while(a->type == typeOpr)
		{
			i++;
			a = a->opr.op[1];
		}
		return i;
	}
	
	int num_out()
	{
		return numOut;
	}



	int add_subfunctions(std::map<std::string, SmartPtr<LUAParserClass> > &subfunctions);
	int add_subfunction(std::string name, std::map<std::string, SmartPtr<LUAParserClass> > &subfunctions);

    int add_subfunctions(std::set<std::string> &knownFunctions, std::stringstream &declarations, std::stringstream &definitions);
    static int addfunctionC(std::string name, std::set<std::string> &knownFunctions, std::stringstream &declarations, std::stringstream &definitions);
    
    void getVar(int i, std::ostream &out);
    int parse_luaFunction(const char *name);
	
    int declare(std::ostream &out);
    int createC_inline(std::ostream &out);
    
    int createRT(nodeType *a, std::ostream &out, const char **rt, int nr, int indent);
    

    int createVMSub(VMAdd &vm, std::map<std::string, SmartPtr<VMAdd> > &subVM);
    int createVMHeader(VMAdd &vm);
    int createVM(nodeType *p, VMAdd &vm,  std::map<std::string, SmartPtr<VMAdd> > &subVM);

    int createVM(VMAdd &vm);

    int	createC(nodeType *p, std::ostream &out, int indent);
    int createJITSG(std::ostream &out, eReturnType r, std::set<std::string> &subfunctions);
	int	createLUA(nodeType *p, std::ostream &out);
	void reduce();
	
	int createC(std::ostream &out);
	int createLUA(std::ostream &out);

	nodeType *reduce(nodeType *p);
	
	/////////
	void print_locals(std::ostream &out);
	void print_globals(std::ostream &out);
	
	/////////////
	static nodeType *con(double value)
	{
		nodeType *p;

		/* allocate node */
		if ((p = new nodeType) == NULL)
			yyerror("out of memory");

		/* copy information */
		p->type = typeCon;
		p->con.value = value;

		return p;
	}
    
    /////////////
	static nodeType *opr0(int i)
	{
		nodeType *p;
		p = new nodeType;
		p->type = typeOpr;
		p->opr.oper = i;
		p->opr.nops = 0;
		return p;
	}

	static nodeType *id(int i)
	{
		nodeType *p;

		/* allocate node */
		if ((p = new nodeType) == NULL)
			yyerror("out of memory");

		/* copy information */
		p->type = typeId;
		p->id.i = i;

		return p;
	}
	
	
	void setRet(int i)
	{
		if(numOut == -1)
			numOut = i;
		else if(numOut != i)
			err << "different return values: was " << numOut << ", now " << i << "\n";
	}

	nodeType *opr1(int oper, nodeType *op)
	{
		nodeType *p;
		p = new nodeType;

		p->type = typeOpr;
		p->opr.oper = oper;
		p->opr.nops = 1;
		p->opr.op = new nodeType*[1];
		p->opr.op[0] = op;
		
		if(oper == 'R')
		{
			nodeType *aa = op;
			int i=1;
			while(aa->type == typeOpr && aa->opr.oper == ',')
			{
				i++;
				aa = aa->opr.op[1];
			}
			setRet(i);			
		}
		return p;
	}
    
	static nodeType *opr2(int oper, nodeType *op1, nodeType *op2)
	{
		nodeType *p;
		p = new nodeType;

		p->type = typeOpr;
		p->opr.oper = oper;
		p->opr.nops = 2;
		p->opr.op = new nodeType*[2];
		p->opr.op[0] = op1;
		p->opr.op[1] = op2;
		return p;
	}
    
    static nodeType *opr3(int oper, nodeType *op1, nodeType *op2, nodeType *op3)
	{
		nodeType *p;
		p = new nodeType;

		p->type = typeOpr;
		p->opr.oper = oper;
		p->opr.nops = 3;
		p->opr.op = new nodeType*[3];
		p->opr.op[0] = op1;
		p->opr.op[1] = op2;
        p->opr.op[2] = op3;
		return p;
	}
	
	static nodeType *forOp(nodeType *_var, nodeType *_start, nodeType *_stop, nodeType *_step, nodeType *_expr)
	{
		nodeType *p;
		p = new nodeType;

		p->type = typeOpr;
		p->opr.oper = LUAPARSER_FOR;
		p->opr.nops = 5;
		p->opr.op = new nodeType*[5];
		p->opr.op[0] = _var;
		p->opr.op[1] = _start;
        p->opr.op[2] = _stop;
		p->opr.op[3] = _step;
		p->opr.op[4] = _expr;
		return p;
	}

    nodeType *function(nodeType *op1, nodeType *op2)
    {
//        std::cout << "local function " << id2variable[op1->id.i] << "\n";
        localFunctions.insert(op1->id.i);
        return opr2('C', op1, op2);
    }

	static void freeNode(nodeType *p)
	{
		int i;

		if (!p) return;
		if (p->type == typeOpr)
		{
			for (i = 0; i < p->opr.nops; i++)
				delete p->opr.op[i];
			delete[] p->opr.op;
		}
		delete p;
	}
	
};

} // namespace ug
#endif	/* LUAParserClass_H */

