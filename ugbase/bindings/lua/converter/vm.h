/*
 * vm.h
 *
 *  Created on: 12.06.2013
 *      Author: mrupp
 */

#ifndef VM_H_
#define VM_H_

#include "common/log.h"
#include <vector>
#include "parser_node.h"

namespace ug{

class VMAdd
{

	std::vector<char> vmBuf;

	enum VMInstruction
	{
		PUSH_CONSTANT=0,
		PUSH_VAR,
		JMP_IF_FALSE,
		JMP,
		OP_UNARY,
		OP_BINARY,
		ASSIGN,
		OP_RETURN
	};

	void serialize(VMInstruction inst)
	{
		int i=inst;
		serialize(i);
	}


	void deserialize(size_t &p, VMInstruction &instr)
	{
		int i = *((int*)&vmBuf[p]);
		instr = (VMInstruction) i;
		p += sizeof(int);
	}

	void serialize(char c)
	{
		vmBuf.push_back(c);
	}

	void serialize(int d)
	{
		char *tmp = (char*)&d;
		for(size_t i=0; i<sizeof(d); i++)
			serialize(tmp[i]);
	}

	void deserialize(size_t &p, int &i)
	{
		i = *((int*)&vmBuf[p]);
		p += sizeof(int);
	}

	void serialize(double d)
	{
		char *tmp = (char*)&d;
		for(size_t i=0; i<sizeof(d); i++)
			serialize(tmp[i]);
	}

	void deserialize(size_t &p, double &d)
	{
		d = *((double*)&vmBuf[p]);
		p += sizeof(double);
	}


	template<typename T>
	void print_unary(const char *desc, size_t &i)
	{
		T t;
		deserialize(i, t);
		UG_LOG(desc << " " << t << "\n");
	}

	void print_op(size_t &i)
	{
		int op;
		deserialize(i, op);
		switch(op)
		{
		case MATH_COS: 		UG_LOG("cos\n"); break;
		case MATH_SIN: 		UG_LOG("sin\n"); break;
		case MATH_EXP: 		UG_LOG("exp\n"); break;
		case MATH_ABS: 		UG_LOG("abs\n"); break;
		case MATH_LOG: 		UG_LOG("log\n"); break;
		case MATH_LOG10: 	UG_LOG("log10\n"); break;
		case MATH_SQRT: 		UG_LOG("sqrt\n"); break;
		case MATH_FLOOR: 		UG_LOG("floor\n"); break;
		case MATH_CEIL: 		UG_LOG("ceil\n"); break;
		case '+': 		UG_LOG("+\n"); break;
		case '-': 		UG_LOG("-\n"); break;
		case '*': 		UG_LOG("*\n"); break;
		case '/': 		UG_LOG("/\n"); break;
		case '<': 		UG_LOG("<\n"); break;
		case '>': 		UG_LOG(">\n"); break;
		case GE: 		UG_LOG("GE\n"); break;
		case LE: 		UG_LOG("LE\n"); break;
		case NE: 		UG_LOG("NE\n"); break;
		case EQ: 		UG_LOG("EQ\n"); break;
		case AND: 		UG_LOG("AND\n"); break;
		case OR: 		UG_LOG("OR\n"); break;
		case MATH_POW: 		UG_LOG("pow\n"); break;
		case MATH_MIN: 		UG_LOG("min\n"); break;
		case MATH_MAX: 		UG_LOG("max\n"); break;
		}
	}

public:
	void push(double constant)
	{
		UG_LOG("POS " << get_pos() << "\n");
		UG_LOG("PUSH_CONSTANT " << constant << "\n");
		serialize((int)PUSH_CONSTANT);
		serialize(constant);
	}

	void push_var(int i)
	{
		UG_LOG("PUSH_VAR " << i << "\n");
		serialize((int)PUSH_VAR);
		serialize(i);
		UG_LOG("POS " << get_pos() << "\n");
	}

	int jmp_if_false()
	{
		UG_LOG("JMP_IF_FALSE\n");
		return jump(JMP_IF_FALSE);
	}

	int jmp()
	{
		UG_LOG("JMP\n");
		return jump(JMP);
	}

	int get_pos()
	{
		return vmBuf.size();
	}

	void unary(int oper)
	{
		UG_LOG("UNARY OP " << oper << "\n");
		serialize((int)OP_UNARY);
		serialize(oper);
		UG_LOG("POS " << get_pos() << "\n");
	}

	void binary(int oper)
	{
		UG_LOG("BINARY OP " << oper << "\n");
		serialize((int)OP_BINARY);
		serialize(oper);
		UG_LOG("POS " << get_pos() << "\n");
	}

	void assign(int v)
	{
		UG_LOG("ASSIGN " << v << "\n");
		serialize((int)ASSIGN);
		serialize(v);
		UG_LOG("POS " << get_pos() << "\n");
	}

	void adjust_jmp_pos(int iPos, int jmpPos)
	{
		UG_LOG("adjusting jmp pos in " << iPos << " to " << jmpPos << "\n");
		int *p = (int *)&vmBuf[iPos];
		*p = jmpPos;
	}

	void ret()
	{
		serialize((int)OP_RETURN);
	}

	void print()
	{
		for(size_t i=0; i<vmBuf.size(); )
		{
			VMInstruction instr;
			deserialize(i, instr);
			UG_LOG(i << "	");
			switch(instr)
			{
				case PUSH_CONSTANT:
					print_unary<double>("PUSH_CONSTANT", i);
					break;

				case PUSH_VAR:
					print_unary<int>("PUSH_VAR", i);
					break;

				case JMP_IF_FALSE:
					print_unary<int>("JMP_IF_FALSE", i);
					break;
				case JMP:
					print_unary<int>("JMP", i);
					break;

				case OP_UNARY:
				case OP_BINARY:
					print_op(i);
					break;
				case ASSIGN:
					print_unary<int>("ASSIGN", i);
					break;
				case OP_RETURN:
					UG_LOG("RETURN\n");
					break;
				default:
					UG_LOG(((int)instr) << "?\n");
			}
		}
	}

	int jump(VMInstruction instr)
	{
		serialize(instr);
		int jmpPos = get_pos();
		serialize(jmpPos);
		UG_LOG("jump pos is " << jmpPos << "\n");
		return jmpPos;
	}

};

}
#endif /* VM_H_ */
