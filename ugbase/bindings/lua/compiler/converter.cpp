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

#include "bridge/bridge.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_stack_check.h"
#include "bindings/lua/info_commands.h"
#include "lua_compiler.h"
#include "lua_parser_class.h"
#include "vm.h"

using namespace std;
namespace ug{
namespace bridge {
	

int convert(const char *functionName)
{
    LUAParserClass parser;
    int ret = parser.parse_luaFunction(functionName);
    if(ret == LUAParserClass::LUAParserError)
	{
		UG_LOG("parsing " << functionName << " failed: reduced LUA parser failed.\n");
		return false;
	}
	if(ret == LUAParserClass::LUAParserIgnore)
	{
		UG_LOG("parsing " << functionName << " : Found --LUACompiler:ignore.\n");
		return false;
	}
    
    if(ret== LUAParserClass::LUAParserOK)
    {
    	UG_LOG("Created C code:\n")
    	parser.createC(cout);
    }

	return 0;	
}

int convertVM(const char *functionName)
{
	LUAParserClass parser;
	int ret = parser.parse_luaFunction(functionName);
	if(ret == LUAParserClass::LUAParserError)
	{
		UG_LOG("parsing " << functionName << " failed: reduced LUA parser failed.\n");
		return false;
	}
	if(ret == LUAParserClass::LUAParserIgnore)
	{
		UG_LOG("parsing " << functionName << " : Found --LUACompiler:ignore.\n");
		return false;
	}

	if(ret== LUAParserClass::LUAParserOK)
	{
		VMAdd vm;
		UG_LOG("Created VM code:\n")
		parser.createVM(vm);
		parser.print_variable_names();
		vm.print();
	}

	return 0;
}


bool RegisterConverter(Registry &reg, const char* parentGroup)
{
	stringstream grpSS; grpSS << parentGroup << "/Converter";
	std::string grp = grpSS.str();

	//try
	{
		reg.add_function("LUA2C_convertC", &convert, grp.c_str());		
		reg.add_function("LUA2C_convertVM", &convertVM, grp.c_str());
	}
	//UG_REGISTRY_CATCH_THROW(grp);

	return true;
}


}}
