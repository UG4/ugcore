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

#ifndef VM_H_
#define VM_H_

#include "common/log.h"
#include <vector>
#include "parser_node.h"
#include "common/assert.h"
#include "parser.hpp"
#include "common/error.h"

namespace ug {

//////////////////////////////////////////
/// --> documentation in vm.doxygen <--///
//////////////////////////////////////////
class VMAdd
{
private:
	std::vector<char> vmBuf;
	std::string m_name;
	std::vector<double> variables;
	size_t m_nrOut, m_nrIn;
	std::vector<SmartPtr<VMAdd> > subfunctions;

	enum VMInstruction
	{
		PUSH_CONSTANT=0,
		PUSH_VAR,
		JMP_IF_FALSE,
		JMP,
		OP_UNARY,
		OP_BINARY,
		ASSIGN,
		OP_RETURN,
		OP_CALL
	};



	void serializeVMInstr(VMInstruction inst)
	{
		int i=inst;
		serializeInt(i);
	}


	inline void deserializeVMInstr(size_t &p, VMInstruction &instr)
	{
		int i = *((int*)&vmBuf[p]);
		instr = (VMInstruction) i;
		p += sizeof(int);
	}

	inline void serializeChar(char c)
	{
		vmBuf.push_back(c);
	}

	inline void serializeInt(int d)
	{
		char *tmp = (char*)&d;
		for(size_t i=0; i<sizeof(d); i++)
			serializeChar(tmp[i]);
	}

	inline void deserializeInt(size_t &p, int &i)
	{
		i = *((int*)&vmBuf[p]);
		p += sizeof(int);
	}

	inline void serializeDouble(double d)
	{
		char *tmp = (char*)&d;
		for(size_t i=0; i<sizeof(d); i++)
			serializeChar(tmp[i]);
	}

	inline void deserializeDouble(size_t &p, double &d)
	{
		d = *((double*)&vmBuf[p]);
		p += sizeof(double);
	}


	void print_unaryDouble(const char *desc, size_t &i)
	{
		double t;
		deserializeDouble(i, t);
		UG_LOG(desc << " " << t << "\n");
	}
	void print_unaryInt(const char *desc, size_t &i)
	{
		int t;
		deserializeInt(i, t);
		UG_LOG(desc << " " << t << "\n");
	}

	void print_op(size_t &i)
	{
		int op;
		deserializeInt(i, op);
		switch(op)
		{
		case LUAPARSER_MATH_COS: 		UG_LOG("cos\n"); break;
		case LUAPARSER_MATH_SIN: 		UG_LOG("sin\n"); break;
		case LUAPARSER_MATH_EXP: 		UG_LOG("exp\n"); break;
		case LUAPARSER_MATH_ABS: 		UG_LOG("abs\n"); break;
		case LUAPARSER_MATH_LOG: 		UG_LOG("log\n"); break;
		case LUAPARSER_MATH_LOG10: 	UG_LOG("log10\n"); break;
		case LUAPARSER_MATH_SQRT: 		UG_LOG("sqrt\n"); break;
		case LUAPARSER_MATH_FLOOR: 		UG_LOG("floor\n"); break;
		case LUAPARSER_MATH_CEIL: 		UG_LOG("ceil\n"); break;
		case '+': 		UG_LOG("+\n"); break;
		case '-': 		UG_LOG("-\n"); break;
		case '*': 		UG_LOG("*\n"); break;
		case '/': 		UG_LOG("/\n"); break;
		case '<': 		UG_LOG("<\n"); break;
		case '>': 		UG_LOG(">\n"); break;
		case LUAPARSER_GE: 		UG_LOG("GE >=\n"); break;
		case LUAPARSER_LE: 		UG_LOG("LE <=\n"); break;
		case LUAPARSER_NE: 		UG_LOG("NE ~=\n"); break;
		case LUAPARSER_EQ: 		UG_LOG("EQ ==\n"); break;
		case LUAPARSER_AND: 		UG_LOG("AND\n"); break;
		case LUAPARSER_OR: 		UG_LOG("OR\n"); break;
		case LUAPARSER_MATH_POW: 		UG_LOG("pow\n"); break;
		case LUAPARSER_MATH_MIN: 		UG_LOG("min\n"); break;
		case LUAPARSER_MATH_MAX: 		UG_LOG("max\n"); break;
		}
	}

public:
	VMAdd()
	{
			m_name = "unknown";
	}
	void set_name(std::string name)
	{
		m_name = name;
	}

	void push(double constant)
	{
//		UG_LOG("POS " << get_pos() << "\n");
//		UG_LOG("PUSH_CONSTANT " << constant << "\n");
		serializeInt(PUSH_CONSTANT);
		serializeDouble(constant);
	}

	void push_var(int i)
	{
//		UG_LOG("PUSH_VAR " << i << "\n");
		serializeInt(PUSH_VAR);
		serializeInt(i);
//		UG_LOG("POS " << get_pos() << "\n");
	}

	int jmp_if_false()
	{
//		UG_LOG("JMP_IF_FALSE\n");
		return jump(JMP_IF_FALSE);
	}

	int jmp()
	{
//		UG_LOG("JMP\n");
		return jump(JMP);
	}

	int get_pos()
	{
		return vmBuf.size();
	}

	void unary(int oper)
	{
//		UG_LOG("UNARY OP " << oper << "\n");
		serializeInt(OP_UNARY);
		serializeInt(oper);
//		UG_LOG("POS " << get_pos() << "\n");
	}

	void binary(int oper)
	{
//		UG_LOG("BINARY OP " << oper << "\n");
		serializeInt(OP_BINARY);
		serializeInt(oper);
//		UG_LOG("POS " << get_pos() << "\n");
	}

	void assign(int v)
	{
//		UG_LOG("ASSIGN " << v << "\n");
		serializeInt(ASSIGN);
		serializeInt(v);
//		UG_LOG("POS " << get_pos() << "\n");
	}

	void call(SmartPtr<VMAdd> subfunction)
	{
		size_t i;
		for(i=0; i<subfunctions.size(); i++)
			if(subfunctions[i] == subfunction)
				break;
		if(i == subfunctions.size())
			subfunctions.push_back(subfunction);
		serializeInt(OP_CALL);
		serializeInt(i);
	}

	void adjust_jmp_pos(int iPos, int jmpPos)
	{
//		UG_LOG("adjusting jmp pos in " << iPos << " to " << jmpPos << "\n");
		int *p = (int *)&vmBuf[iPos];
		*p = jmpPos;
	}

	void ret()
	{
		serializeInt(OP_RETURN);
	}

	void print_short()
	{
		UG_LOG("function " << m_name << ", " << m_nrIn << " inputs, " << m_nrOut <<
				" outputs, "<< variables.size() << " variables, " << subfunctions.size() << " subfunctions");
	}

	void print()
	{
		print_rec(0);
	}

	void print_rec(int level)
	{
		if(level > 5) { UG_LOG("\n... aborting recursion (potential infinite loop)\n"); }
		print_short();
		UG_LOG("\n");
		for(size_t i=0; i<vmBuf.size(); )
		{
			UG_ASSERT(i<vmBuf.size(), i);
			UG_LOG(i << "	");
			VMInstruction instr;
			deserializeVMInstr(i, instr);

			switch(instr)
			{
				case PUSH_CONSTANT:
					print_unaryDouble("PUSH_CONSTANT", i);
					break;

				case PUSH_VAR:
					print_unaryInt("PUSH_VAR", i);
					break;

				case JMP_IF_FALSE:
					print_unaryInt("JMP_IF_FALSE", i);
					break;
				case JMP:
					print_unaryInt("JMP", i);
					break;

				case OP_UNARY:
					UG_LOG("OP_UNARY "); print_op(i);
					break;
				case OP_BINARY:
					UG_LOG("OP_BINARY "); print_op(i);
					break;
				case ASSIGN:
					print_unaryInt("ASSIGN", i);
					break;
				case OP_RETURN:
					UG_LOG("RETURN\n");
					break;

				case OP_CALL:
				{
					int varI;
					deserializeInt(i, varI);

					UG_COND_THROW(varI >= (int)subfunctions.size(), i);
					UG_LOG("CALL to subfunction " << varI << ": ");
					subfunctions[varI]->print_short();
					UG_LOG("\n");
					break;
				}

				default:
					UG_LOG(((int)instr) << "?\n");
			}
		}
		if(subfunctions.size() > 0)
		{
			UG_LOG("<< SUBFUNCTIONS OF " << m_name << ":\n");
			for(size_t i=0; i<subfunctions.size(); i++)
			{
				UG_LOG("---- SUBFUNCTION " << i << " of " << m_name << " ----\n");
				subfunctions[i]->print_rec(level+1);
			}
			UG_LOG("   SUBFUNCTIONS OF " << m_name << " >>\n");
		}
	}

	int jump(VMInstruction instr)
	{
		serializeVMInstr(instr);
		int jmpPos = get_pos();
		serializeInt(jmpPos);
//		UG_LOG("jump pos is " << jmpPos << "\n");
		return jmpPos;
	}



	void set_in_out(size_t nrIn,size_t nrOut)
	{
		m_nrOut = nrOut;
		m_nrIn = nrIn;
	}
	void set_nr_of_variables(size_t nr)
	{
		variables.resize(nr);
	}

	inline void execute_unary(size_t &i, double &v)
	{
		int op;
//		UG_LOG("unary op " << v << "\n");
		deserializeInt(i, op);
		switch(op)
		{
			case LUAPARSER_MATH_COS: v = cos(v);	break;
			case LUAPARSER_MATH_SIN: v = sin(v);	break;
			case LUAPARSER_MATH_EXP: v = exp(v);	break;
			case LUAPARSER_MATH_ABS: v = fabs(v); 	break;
			case LUAPARSER_MATH_LOG: v = log(v); 	break;
			case LUAPARSER_MATH_LOG10: v = log10(v); break;
			case LUAPARSER_MATH_SQRT:  v = sqrt(v); 	break;
			case LUAPARSER_MATH_FLOOR: v = floor(v); break;
			case LUAPARSER_MATH_CEIL: v = ceil(v); break;
		}
	}

	inline void execute_binary(size_t &i, double *stack, int SP)
	{
		double &a = stack[SP-2];
		double &b = stack[SP-1];
//		UG_LOG("binary op " << a << " op " << b << "\n");

		int op;
		deserializeInt(i, op);
		switch(op)
		{
			case '+': 	a = b+a;	break;
			case '-': 	a = b-a;	break;
			case '*': 	a = b*a;	break;
			case '/': 	a = b/a;	break;
			case '<': 	a = (b < a) ? 1.0 : 0.0; break;
			case '>': 	a = (b > a) ? 1.0 : 0.0; break;
			case LUAPARSER_GE: 	a = (b >= a) ? 1.0 : 0.0; break;
			case LUAPARSER_LE: 	a = (b <= a) ? 1.0 : 0.0; break;
			case LUAPARSER_NE: 	a = (b != a) ? 1.0 : 0.0; break;
			case LUAPARSER_EQ: 	a = (b == a) ? 1.0 : 0.0; break;
			case LUAPARSER_AND: 	a = (a != 0.0 && b != 0.0) ? 1.0 : 0.0; break;
			case LUAPARSER_OR: 	a = (a != 0 || b != 0) ? 1.0 : 0.0; break;
			case LUAPARSER_MATH_POW: 	a = pow(b, a); break;
			case LUAPARSER_MATH_MIN: 	a = (b < a) ? a : b; break;
			case LUAPARSER_MATH_MAX: 	a = (b > a) ? a : b; break;
		}
	}

	double call(double *stack, int &SP)
	{
//		for(int j=0; j<variables.size(); j++)
//		{UG_LOG("var[[" << j << "] = " << variables[j] << "\n");}

		double varD;
		int varI;
		size_t i=0;
		VMInstruction instr;
		while(true)
		{
//			UG_LOG("IP =  " << i << ", SP = " << SP);
			deserializeVMInstr(i, instr);
//			UG_LOG("OP =  " << (int)instr << "\n");
//			for(int j=0; j<SP; j++)
//				{UG_LOG("SP[" << j << "] = " << stack[j] << ", ");}
//			UG_LOG("\n");

			switch(instr)
			{
				case PUSH_CONSTANT:

					deserializeDouble(i, varD);
//					UG_LOG("PUSH CONSTANT " << varD << "\n");
					stack[SP++] = varD;
					break;

				case PUSH_VAR:
					deserializeInt(i, varI);
					stack[SP++] = variables[varI-1];
//					UG_LOG("PUSH VAR " << varI << "\n");
					break;

				case JMP_IF_FALSE:
					deserializeInt(i, varI);
//					UG_LOG("JMP IF FALSE " << varI << "\n");
					SP--;
					if(stack[SP] == 0.0)
						i = varI;
					break;
				case JMP:
					deserializeInt(i, varI);
//					UG_LOG("JMP " << varI << "\n");
					i = varI;
					break;

				case OP_UNARY:
					UG_ASSERT(SP>0, SP);
					execute_unary(i, stack[SP-1]);
					break;

				case OP_BINARY:
					UG_ASSERT(SP>1, SP);
					execute_binary(i, stack, SP);
					SP--;
					break;

				case ASSIGN:
					deserializeInt(i, varI);
					SP--;
//					UG_LOG("ASSIGN " << varI << " = " << stack[SP] << "\n");
					variables[varI-1] = stack[SP];
					break;

				case OP_RETURN:
//					UG_LOG("RETURN\n");
					UG_ASSERT(SP == (int)m_nrOut, "stack pointer is not nrOut =" << m_nrOut << ", instead " << SP << " ?")
					return stack[0];
					break;

				case OP_CALL:
				{
					deserializeInt(i, varI);
//					UG_LOG("call IP =  " << i << ", SP = " << SP << "\n");
					SmartPtr<VMAdd> sub = subfunctions[varI];
					sub->call_sub(stack, SP);
//					UG_LOG("call IP =  " << i << ", SP = " << SP << "\n");
					break;
				}
				default:
					UG_ASSERT(0, "IP: " << i << " op " << ((int)instr) << " ?\n");
			}
		}
	}

	double call_sub(double *stack, int &SP)
	{
		for(size_t i=0; i<m_nrIn; i++)
			variables[i] = stack[SP+i-1];
		SP -= m_nrIn;
		return call(stack, SP);
	}

	int execute(double *ret, const double *in)
	{
		double stack[255];
		int SP=0;

		for(size_t i=0;i<m_nrIn; i++)
			variables[i] = in[i];
		call(stack, SP);
		UG_ASSERT(SP == (int)m_nrOut, SP << " != " << m_nrOut);
		for(size_t i=0; i<m_nrOut; i++)
			ret[i] = stack[i];

		return 1;
	}

	double call()
	{
		double stack[255];
		int SP=0;
		return call(stack, SP);
	}

	double call(double a, double b, double c)
	{
		UG_ASSERT(m_nrOut == 1, m_nrOut);
		UG_ASSERT(m_nrIn == 3, m_nrIn);
		variables[0] = a;
		variables[1] = b;
		variables[2] = c;
		return call();
	}

	double call(double a, double b)
	{
		UG_ASSERT(m_nrOut == 1, m_nrOut);
		UG_ASSERT(m_nrIn == 2, m_nrIn);
		variables[0] = a;
		variables[1] = b;
		return call();
	}

	double call(double a)
	{
		UG_ASSERT(m_nrOut == 1, m_nrOut);
		UG_ASSERT(m_nrIn == 1, m_nrIn);
		variables[0] = a;
		return call();
	}

	size_t num_out()
	{
		return m_nrOut;
	}
	size_t num_in()
	{
		return m_nrIn;
	}
};

}
#endif