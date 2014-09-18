/*
 * \file	converter.cpp
 * \author	Martin Rupp
 *
 * Created on 20. November 2012, 10:16
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
