/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
#include "bridge/bridge.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_stack_check.h"
#include "bindings/lua/info_commands.h"
#include "common/util/string_util.h"
#include "lua_compiler_debug.h"
#include "common/profiler/profiler.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

using namespace std;

namespace ug {
    

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

int LUAParserClass::get_id_for_name(const char*name)
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

void LUAParserClass::print_variable_names()
{
	for(std::map<size_t, std::string>::iterator it = id2variable.begin(); it != id2variable.end(); ++it)
	{
		UG_LOG((*it).first << ": " << (*it).second);
	}
}

void LUAParserClass::getVar(int i, ostream &out)
{
	PROFILE_BEGIN_GROUP(LUAParserClass_getVar, "LUA2C");
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

int LUAParserClass::createRT(nodeType *a, ostream &out, const char **rt, int nr, int indent)
{
	PROFILE_BEGIN_GROUP(LUAParserClass_createRT, "LUA2C");
    int i=0;
    while(a->type == typeOpr && a->opr.oper == ',')
    {
    	out << repeat('\t', indent);
        out << rt[i++] << " = ";
        createC(a->opr.op[0], out, indent);
        out << ";\n";
        a = a->opr.op[1];
    }
    out << repeat('\t', indent);
    out << rt[i++] << " = ";
    createC(a, out, indent);
    out << ";\n";
    out << repeat('\t', indent);
    return true;
}


void LUAParserClass::reduce()
{
	PROFILE_BEGIN_GROUP(LUAParserClass_reduce, "LUA2C");
	PROFILE_FUNC();
	for(size_t i=0; i<nodes.size(); i++)
		nodes[i] = reduce(nodes[i]);
}

#ifdef UG_PARALLEL
string GetFileLinesParallel(string filename, size_t fromline, size_t toline, bool includeLineNumbers, const pcl::ProcessCommunicator &pc = pcl::ProcessCommunicator())
{
	PCL_PROFILE_FUNC();
	string lines;
	BinaryBuffer buf;
	if(pc.is_proc_id(0))
	{
		lines = GetFileLines(filename.c_str(), fromline, toline, includeLineNumbers);
		Serialize(buf, lines);
	}

	pc.broadcast(buf);

	if(!pc.is_proc_id(0))
		Deserialize(buf, lines);

	return lines;
}
#endif

int LUAParserClass::parse_luaFunction(LuaFunctionHandle handle)
{
    lua_State* L = script::GetDefaultLuaState();
	lua_rawgeti(L, LUA_REGISTRYINDEX, handle.ref);

	if(!lua_isfunction(L, -1))
	{
		UG_DLOG(DID_LUACOMPILER, 1, "LUA Script function " << "__unknown__lua__function__by__handle" << " not found\n");
		lua_pop(L, 1);
		return false;
	}

	return parse_luaFunction_StackTop("__unknown__lua__function__by__handle");
}

int LUAParserClass::parse_luaFunction(const char *functionName)
{
	PROFILE_BEGIN_GROUP(LUAParserClass_parse_luaFunction, "LUA2C");
	PROFILE_FUNC();
    lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);

	//UG_LOG("LUA2C::create(" << functionName << ")\n");

	bridge::GetLuaNamespace(L, functionName);

	if(!lua_isfunction(L, -1))
	{
		UG_DLOG(DID_LUACOMPILER, 1, "LUA Script function " << functionName << " not found\n");
		lua_pop(L, 1);
		return false;			
	}

	return parse_luaFunction_StackTop(functionName);
}

int LUAParserClass::parse_luaFunction_StackTop(const char *functionName)
{
    lua_State* L = script::GetDefaultLuaState();

    lua_pushvalue(L, -1);
	lua_Debug ar;
	lua_getinfo(L, ">Snlu", &ar);
	if(!ar.source)
	{
		UG_DLOG(DID_LUACOMPILER, 1, "no source found\n")
		lua_pop(L, 1);
		//UG_ASSERT(0, "no source found");
		return false;			
	}
	const char *src = ar.source[0]=='@' ? ar.source+1 : ar.source;

#ifdef UG_PARALLEL
	string str = GetFileLinesParallel(src, ar.linedefined, ar.lastlinedefined, false);
#else
	string str = GetFileLines(src, ar.linedefined, ar.lastlinedefined, false);
#endif

	lua_pop(L, 1);

	// check for --LUACompiler:ignore
	if(ar.linedefined > 0)
	{
		string str2 = GetFileLines(src, ar.linedefined-1, ar.linedefined, false);
		if(str2.find("--lua2c:ignore") != std::string::npos || str2.find("--LUACompiler:ignore") != std::string::npos)
			return LUAParserIgnore;
	}
	if(str.find("--lua2c:ignore") != std::string::npos || str.find("--LUACompiler:ignore") != std::string::npos)
		return LUAParserIgnore;

	//UG_DLOG("The function:\n"<<str<<"\n");


    iLineAdd = ar.linedefined;
    filename = src;
    filename = FilenameWithoutPath(filename);
    
	parse(str.c_str());
    
    if(has_errors())
	{
    	UG_LOG("-----[ LUACompiler parsing error for function " << functionName << ":  -----\n");
    	UG_LOG("-- by adding --LUACompiler:ignore to " << functionName << ", LUACompiler will ignore this function --\n");
    	UG_LOG(src << " " << ar.linedefined << " - " << ar.lastlinedefined << " : \n");
    	UG_LOG(GetFileLines(src, ar.linedefined, ar.lastlinedefined, true));
    	UG_LOG("\n----- Parsing errors:\n");
    	UG_LOG(err.str());
    	UG_LOG("-----]\n");

		//UG_ASSERT(0, parser.err.str());
		return LUAParserError;
	}
    return LUAParserOK;
}

int LUAParserClass::add_subfunctions(set<string> &knownFunctions, stringstream &declarations, stringstream &definitions)
{
    //UG_LOG("nr of subfunctions: " << localFunctions.size() << ".\n");        
    for(set<size_t>::iterator it = localFunctions.begin(); it != localFunctions.end(); ++it)
    {
    	int ret = addfunctionC(id2variable[*it], knownFunctions, declarations, definitions);
    	if(ret != LUAParserOK)
    		return ret;
    }
    return LUAParserOK;
}





void LUAParserClass::print_locals(ostream &out)
{
	out << "local variables:\n";
	for(map<string, size_t>::iterator it = variables.begin(); it != variables.end(); ++it)
		if(is_local((*it).second))
			out << "\t" << (*it).first << "\n";

}
void LUAParserClass::print_globals(ostream &out)
{
	out << "global references:\n";
	for(map<string, size_t>::iterator it = variables.begin(); it != variables.end(); ++it)
		if(!is_local((*it).second))
			out << "\t" << (*it).first << "\n";
}

}
