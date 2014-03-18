/*
 * \file	lua2c.cpp
 * \author	Martin Rupp
 *
 * Created on 20. November 2012, 10:16
 */

#include "lua_parser_class.h"
#include "bridge/bridge.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_stack_check.h"
#include "bindings/lua/info_commands.h"
#include "common/util/file_util.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "lua2c.h"
using namespace std;

namespace ug{

extern bool useLua2VM;
extern bool useLua2C;

namespace bridge {
    
DebugID DID_LUA2C("LUA2C");
bool LUA2C::create(const char *functionName)
{
	if(useLua2VM)
		return createVM(functionName);
	else
		return createC(functionName);
}

bool LUA2C::createC(const char *functionName)
{
	UG_LOG("parsing " << functionName << "... ");
	try{
		m_f=NULL;
		LUAParserClass parser;
		if(parser.parse_luaFunction(functionName) == false)
			return 0;
		//parser.reduce();
	
		string p = PathProvider::get_path(ROOT_PATH) + "/bin/lua2c_tmp/";

		if(!DirectoryExists(p))
			CreateDirectory(p);

		fstream out((p + "lua2c_output.c").c_str(), fstream::out);
	
		out << "#include <math.h>\n";
		out << "#define true 1\n";
		out << "#define false 0\n";
		parser.createC(out);

		m_iIn = parser.num_in();
		m_iOut = parser.num_out();
		out.close();

		string c1s=string("gcc -fpic -O3 -c ") + p + "lua2c_output.c -o " + p + "lua2c_output.o";
		UG_DLOG(DID_LUA2C, 2, "compiling line: " << c1s << "\n");
		if(system(c1s.c_str()) != 0)
		{
			UG_LOG("Error when compiling " << functionName << "\n");
			UG_LOG("compiling line: " << c1s << "\n");
			UG_LOG("--[LUA2C]-----------------------------------------------\n");
			UG_LOG("created C function from LUA function " << functionName << ":\n");
			UG_LOG(GetFileLines((p+"lua2c_output.c").c_str(), 1, -1, true) << "\n");
			UG_LOG("--[LUA2C]-----------------------------------------------\n");
			return false;
		}
	
		bool bTmpFileSuccess=false;
		m_name = functionName;

		m_pDyn = MakeTmpFile(p+string(functionName), ".dylib", bTmpFileSuccess);

		string c2s=string("gcc -dynamiclib ") + p+"lua2c_output.o -o " + m_pDyn.c_str();
		UG_DLOG(DID_LUA2C, 2, "linking line: " << c2s << "\n");
		if(system(c2s.c_str()) != 0)
		{
			UG_LOG("Error when linking " << functionName << "\n");
			UG_LOG("linking line: " << c2s << "\n");
			UG_LOG("--[LUA2C]-----------------------------------------------\n");
			UG_LOG("created C function from LUA function " << functionName << ":\n");
			UG_LOG(GetFileLines((p+"lua2c_output.c").c_str(), 1, -1, true) << "\n");
			UG_LOG("--[LUA2C]-----------------------------------------------\n");
			return false;
		}
		try{
		m_libHandle = OpenLibrary(m_pDyn.c_str());
		}
		catch(std::string error)
		{
			UG_LOG("Error when opening library for function " << functionName << "\n");
			UG_LOG("Error is " << error << "\n");
			return false;
		}
		m_f = (LUA2C_Function) GetLibraryProcedure(m_libHandle, functionName);

		if(m_f !=NULL) { UG_LOG("OK\n"); }
		else { UG_LOG("FAILED.\n"); }
		if(m_f !=NULL)
			bInitialized = true;
		return m_f != NULL;
	}
	catch(...)
	{
		UG_LOG("exception thrown in LUA2C::create(" << functionName << ")\n");
		return false;
	}
	
}

bool LUA2C::createVM(const char *functionName)
{
	m_name = functionName;

	LUAParserClass parser;
	try{
		if(parser.parse_luaFunction(functionName) == false)
		{
			UG_LOG("parsing " << functionName << "... ");
			UG_LOG("failed.\n");
			return false;
		}

		if(parser.createVM(vm) == false)
		{
			UG_LOG("parsing " << functionName << "... ");
			UG_LOG("failed.\n");
			return false;
		}
	}
	catch(UGError e)
	{
		UG_LOG("parsing " << functionName << "... ");
		UG_LOG("failed:\n" << e.get_stacktrace() << "\n");
		return false;
	}
	catch(...)
	{
		UG_LOG("parsing " << functionName << "... ");
		UG_LOG("failed: Exception.\n");
		return false;
	}
	//UG_LOG(" ok.\n");
//	vm.print();
	m_iIn = vm.num_in();
	m_iOut = vm.num_out();
	bInitialized = true;
	bVM = true;
	return true;
}


LUA2C::~LUA2C()
{
    UG_DLOG(DID_LUA2C, 2, "removing " << m_name << "\n");		
	if(m_libHandle)
	   	CloseLibrary(m_libHandle);
        
    if(m_pDyn.size() > 0)
    {
		string s = string("rm ") + m_pDyn;
		UG_DLOG(DID_LUA2C, 2, s << "\n");
		system(s.c_str());
	}
}


}
}
