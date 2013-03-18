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
#include "lua2c.h"
#include "pclass.h"

using namespace std;
namespace ug{
namespace bridge {
	

int convert(const char *functionName)
{
    pclass parser;
    if(parser.parse_luaFunction(functionName) == false) return 0;
    
    parser.createC(cout);

	return 0;	
}


bool RegisterConverter(Registry &reg, const char* parentGroup)
{
	stringstream grpSS; grpSS << parentGroup << "/Converter";
	std::string grp = grpSS.str();

	//try
	{
		reg.add_function("LUA2C_convertC", &convert, grp.c_str());		
	}
	//UG_REGISTRY_CATCH_THROW(grp);

	return true;
}


}}
