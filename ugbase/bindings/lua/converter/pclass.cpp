/* 
 * File:   pclass.h
 * Author: mrupp
 *
 * Created on 20. November 2012, 10:16
 */

#include "pclass.h"
#include <cassert>
#include <vector>
#include <sstream>
#include "bridge/bridge.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_stack_check.h"
#include "bindings/lua/info_commands.h"


using namespace std;

namespace ug{
    

vector<nodeType*> CommaListToVector(nodeType *a)
{
    vector<nodeType *> v;
    while(a->type == typeOpr)
	{
        v.push_back(a->opr.op[0]);
		a = a->opr.op[1];
	}
    v.push_back(a);
    return v;
}

int pclass::get_id_for_name(const char*name)
{

	size_t s = variables.size();
	size_t &v = variables[string(name)];
	if (s != variables.size())
	{
		v = s;
		id2variable[v] = name;
		//printf("allocated variable %s, idx %d\n", name, v);
	}
	return v;
}

void repeat(ostream &out, const char *c, int times)
{
	for(int i=0; i<times; i++)
		out << c;
}


void pclass::getVar(int i, ostream &out)
{
    if(is_local(i) || is_arg(i))
        out << id2variable[i].c_str();
    else
    {
        if(id2variable[i].compare("true")==0)
            out << 1;
        else if(id2variable[i].compare("false")==0)
            out << 0;
        else
        {
            lua_State* L = ug::script::GetDefaultLuaState();
            out << ug::bridge::LuaGetNumber(L, id2variable[i].c_str(), 0);
        }
    }
}

int pclass::createRT(nodeType *a, ostream &out, const char **rt, int nr, int indent)
{
    int i=0;
    while(a->type == typeOpr && a->opr.oper == ',')
    {
        repeat(out, "\t", indent);
        out << rt[i++] << " = ";
        createC(a->opr.op[0], out, indent);
        out << ";\n";
        a = a->opr.op[1];
    }
    repeat(out, "\t", indent);
    out << rt[i++] << " = ";
    createC(a, out, indent);
    out << ";\n";
    repeat(out, "\t", indent);
    return true;
}
int pclass::createC(nodeType *p, ostream &out, int indent)
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
				case IF:
                {
					repeat(out, "\t", indent);			out << "if(";
					createC(p->opr.op[0], out, 0);
					out << ")\n";
					repeat(out, "\t", indent); 			out << "{\n";
					createC(p->opr.op[1], out, indent+1);
					repeat(out, "\t", indent); 			out << "}\n";
                 
                    nodeType *a = p->opr.op[2];
                    while(a != NULL && a->opr.oper == ELSEIF)
                    {
                        repeat(out, "\t", indent); 			out << "else if(";
                        createC(a->opr.op[0], out, 0);
                        out << ")\n";
                        repeat(out, "\t", indent); 			out << "{\n";
                        createC(a->opr.op[1], out, indent+1);
                        repeat(out, "\t", indent); 			out << "}\n";
                        a = a->opr.op[2];
                    }
                    if(a != NULL)
                    {
                        UG_ASSERT(a->opr.oper == ELSE, a->opr.oper);
                        repeat(out, "\t", indent); 			out << "else\n";
                        repeat(out, "\t", indent); 			out << "{\n";
                        createC(a->opr.op[0], out, indent+1);
                        repeat(out, "\t", indent); 			out << "}\n";
                    }
                                
                    break;
                }
                
				case '=':
                    UG_ASSERT(is_local(p->opr.op[0]->id.i), "global variable " << id2variable[p->opr.op[0]->id.i] << " is read-only");					
					repeat(out, "\t", indent);
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
					out << ")\n";     
                    break;
                }
                    
				case 'R':
				{
					if(returnType == RT_SUBFUNCTION)
                    {
                        nodeType *a = p->opr.op[0];
                        if(a->type == typeOpr && a->opr.oper == ',')
                        {
                            UG_ASSERT(0, "subfunctions may not return more then one value");
                            return false;                        
                        }
                        repeat(out, "\t", indent);
                        out << "return ";
                        createC(a, out, indent);
                        out << ";\n";
                    }
                    else if(returnType == RT_DIFFUSION)
                    {
                        const char *rt[] = {"A11", "A12", "A21", "A22"};
                        createRT(p->opr.op[0], out, rt, 4, indent);
                        
                        out << "\n";
                        repeat(out, "\t", indent); 
                        out << "goto diffusionReturn;\n";
                    }
                    else if(returnType == RT_VELOCITY)
                    {
                        const char *rt[] = {"vx", "vy"};
                        createRT(p->opr.op[0], out, rt, 2, indent);
                        
                        out << "\n";
                        repeat(out, "\t", indent); 
                        out << "goto velocityReturn;\n";
                    }
                    else if(returnType == RT_DIRICHLET || returnType == RT_SOURCE)
                    {
                        const char *rt[] = {"f"};
                        createRT(p->opr.op[0], out, rt, 1, indent);
                        
                        out << "\n";
                        repeat(out, "\t", indent); 
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
                            repeat(out, "\t", indent);
                            out << "LUA2C_ret[" << i++ << "] = ";
                            createC(a->opr.op[0], out, indent);
                            out << ";\n";
                            a = a->opr.op[1];
                        }
                        repeat(out, "\t", indent);
                        out << "LUA2C_ret[" << i++ << "] = ";
                        createC(a, out, indent);
                        out << ";\n";
                        repeat(out, "\t", indent);
                        out << "return 1;\n";
                    }
                    break;
				}
				
				case TK_FOR:
					repeat(out, "\t", indent);
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
					repeat(out, "\t", indent);
					out << "{\n";
						createC(p->opr.op[4], out, indent+1);
					repeat(out, "\t", indent);
					out << "}\n";
					break;
					
				case TK_BREAK:
					repeat(out, "\t", indent);
					out << "break;\n";
					break;
					
				case UMINUS:
					out << "-(";
					createC(p->opr.op[0], out, 0);
					out << ")";
					break;
                    
                case MATH_PI:
                    out << " MATH_PI ";
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
                        case MATH_COS: out << "cos("; break;
                        case MATH_SIN: out << "sin("; break;
                        case MATH_EXP: out << "exp("; break;
                        case MATH_ABS: out << "fabs("; break;
                        case MATH_LOG: out << "log("; break;
                        case MATH_LOG10: out << "log10("; break;
                        case MATH_SQRT: out << "sqrt("; break;
                        case MATH_FLOOR: out << "floor("; break;
                        case MATH_CEIL: out << "ceil("; break;
                    }                        
					createC(p->opr.op[0], out, 0);
					out << ")";
					break;

                case MATH_POW:
                case MATH_MIN:
                case MATH_MAX:
                    switch (p->opr.oper)
                    {
                        case MATH_POW: out << "pow("; break;
                        case MATH_MIN: out << "min("; break;
                        case MATH_MAX: out << "max("; break;
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
						case GE: out << " >= ";
							break;
						case LE: out << " <= ";
							break;
						case NE: out << " != ";
							break;
						case EQ: out << " == ";
							break;
						case AND: out << " && ";
							break;
						case OR: out << " || ";
							break;
					}
					out << "(";
					createC(p->opr.op[1], out, 0);
					out << ")";
			}
	}
	return 0;
}

int pclass::createLUA(nodeType *p, ostream &out)
{
	if (!p) return 0;
    nodeType *a;
    int i;
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
					i=0;
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
					i=0;
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
void pclass::reduce()
{
	for(size_t i=0; i<nodes.size(); i++)
		nodes[i] = reduce(nodes[i]);
}

int pclass::parse_luaFunction(const char *functionName)
{
    lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);

	//UG_LOG("LUA2C::create(" << functionName << ")\n");

	bridge::GetLuaNamespace(L, functionName);

	if(!lua_isfunction(L, -1))
	{
		UG_LOG("LUA Script function " << functionName << " not found\n");
		lua_pop(L, 1);
		UG_ASSERT(0, "LUA Script function " << functionName << " not found\n");
		return false;			
	}


	lua_pushvalue(L, -1);
	lua_Debug ar;
	lua_getinfo(L, ">Snlu", &ar);
	if(!ar.source)
	{
		UG_LOG("no source found\n")
		lua_pop(L, 1);
		//UG_ASSERT(0, "no source found");
		return false;			
	}
	const char *src = ar.source[0]=='@' ? ar.source+1 : ar.source;

	string str = GetFileLines(src, ar.linedefined, ar.lastlinedefined, false);
	//UG_DLOG("The function:\n"<<str<<"\n");

    lua_pop(L, 1);
    
	parse(str.c_str());
    
    if(has_errors())
	{
		UG_LOG("-- Parsing function " << functionName << ":\n");
		UG_LOG(GetFileLines(src, ar.linedefined, ar.lastlinedefined, true));
		UG_LOG("\nReturned errors:\n");
		UG_LOG(err.str() << "\n");
		//UG_ASSERT(0, parser.err.str());
		return false;
	}
    return true;
}

int pclass::add_subfunctions(set<string> &knownFunctions, stringstream &declarations, stringstream &definitions)
{
    //UG_LOG("nr of subfunctions: " << localFunctions.size() << ".\n");        
    for(set<size_t>::iterator it = localFunctions.begin(); it != localFunctions.end(); ++it)
        if(addfunctionC(id2variable[*it], knownFunctions, declarations, definitions) == false) return false;
    return true;
}

int pclass::addfunctionC(string name, set<string> &knownFunctions, stringstream &declarations, stringstream &definitions)
{
    //UG_LOG("adding " << name << "\n");
    if(knownFunctions.find(name) != knownFunctions.end()) return true;
    knownFunctions.insert(name);
    
    pclass parser;
    if(parser.parse_luaFunction(name.c_str()) == false)
        return false;
    
    if(parser.num_out() != 1)
    {
        UG_LOG("error in subfunction " << name << ": subfunction must have exactly one return value (not " << parser.num_out() << ")\n");
        return false;
    }
    
    parser.returnType = RT_SUBFUNCTION;
    
    parser.declare(declarations); declarations << ";\n";
    
    parser.createC_inline(definitions);
    
    parser.add_subfunctions(knownFunctions, declarations, definitions);
    
    return true;    
}
int pclass::createJITSG(ostream &out, eReturnType r, set<string> &subfunctions)
{	
    int numRet=0;
    switch(r)
    {
        case RT_DIRICHLET: numRet = 1; break;
        case RT_SOURCE: numRet = 1; break;
        case RT_DIFFUSION: numRet = 4; break;
        case RT_VELOCITY: numRet = 2; break;
        default: UG_ASSERT(0,"");
    }
    if(num_out() != numRet)
    {
        UG_LOG("number of return values must be " << numRet << ", not " << num_out() << "\n");
        return false;
    }
    out << "// INSERTED CODE:";	    
    out << "\n// JITSG Generated Code from LUA function " << name << "\n";	
    
    
    for(set<size_t>::iterator it = localFunctions.begin(); it != localFunctions.end(); ++it)
        subfunctions.insert(id2variable[*it]);
    
    // escape variables
    for(map<size_t, string>::iterator it = id2variable.begin(); it != id2variable.end(); ++it)
		if(is_local((*it).first))
			(*it).second = "_" + (*it).second;
    // TODO: also escape functions
    
    // rename arguments to x, y, t
    nodeType *a = args;
    id2variable[a->opr.op[0]->id.i] = "x";
    a = a->opr.op[1];
    id2variable[a->opr.op[0]->id.i] = "y";
    id2variable[a->opr.op[1]->id.i] = "t";
    
    //------ local variables --------
	out << "\t// local variables:\n";
	for(map<size_t, string>::iterator it = id2variable.begin(); it != id2variable.end(); ++it)
		if(is_local((*it).first))
			out << "\tdouble " << (*it).second << ";\n";
	
    returnType = r;
			
	out << "\t// code:\n";
	for(size_t i=0; i<nodes.size(); i++)
		createC(nodes[i], out, 1);
	out << "\n// END JITSG Generated Code from LUA function " << name << "\n";
	switch(r)
	{
		case RT_DIRICHLET: out << "dirichletReturn:\n"; break;
        case RT_SOURCE: out << "sourceReturn:\n"; break;
        case RT_DIFFUSION: out << "diffusionReturn:\n"; break;
        case RT_VELOCITY: out << "velocityReturn:\n"; break;
		default: UG_ASSERT(0,"?");
	}    
	out << ";\n";
	return true;
}

int pclass::createC(ostream &out)
{
     // local functions    
    
    set<string> knownFunctions;
    stringstream declarations;
    stringstream definitions;
    if(add_subfunctions(knownFunctions, declarations, definitions) == false)
    {
        UG_LOG("add_subfunctions failed.\n");        
        return false;
    }
     
    out << "#define MATH_PI 3.1415926535897932384626433832795028841971693\n";
    out << "// inline function declarations\n";
    out << declarations.str() << "\n";
    
    out << "// inline function definitions\n";
    out << definitions.str() << "\n";
    
    // the function
	out << "int " << name << "(";
	out << "double *LUA2C_ret, "; //[" << numOut << "], ";
	out << "double *LUA2C_in)\n"; //[" << numIn << "])\n";
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
	return 0;
}

int pclass::declare(ostream &out)
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

int pclass::createC_inline(ostream &out)
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

int pclass::createLUA(ostream &out)
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

nodeType *
pclass::reduce(nodeType *p)
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
				case IF:
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
				case UMINUS:
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

				case MATH_COS:
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

				case MATH_SIN:
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

				case MATH_EXP:
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
							case GE: return con(a > b);
							case LE: return con(a < b);
							case NE: return con(a != b);
							case EQ: return con(a == b);
							case AND: return con(a != 0.0 && b != 0.0 ? 1.0: 0.0);
							case OR: return con(a != 0.0 || b != 0.0 ? 1.0: 0.0);
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


void pclass::print_locals(ostream &out)
{
	out << "local variables:\n";
	for(map<string, size_t>::iterator it = variables.begin(); it != variables.end(); ++it)
		if(is_local((*it).second))
			out << "\t" << (*it).first << "\n";

}
void pclass::print_globals(ostream &out)
{
	out << "global references:\n";
	for(map<string, size_t>::iterator it = variables.begin(); it != variables.end(); ++it)
		if(!is_local((*it).second))
			out << "\t" << (*it).first << "\n";
}

}
