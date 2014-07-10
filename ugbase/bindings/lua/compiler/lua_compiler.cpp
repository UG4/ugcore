/*
 * \file	LUACompiler.cpp
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
#include "lua_compiler.h"
#include "lua_compiler_debug.h"
#include "common/profiler/profiler.h"
using namespace std;

namespace ug{

extern bool useLua2VM;
extern bool useLuaCompiler;

DebugID DID_LUACOMPILER("LUACompiler");

namespace bridge {
    

bool LUACompiler::create(const char *functionName)
{
	if(useLua2VM)
		return createVM(functionName);
	else
		return createC(functionName);
}

bool LUACompiler::createC(const char *functionName)
{
#ifdef USE_LUA2C
	PROFILE_BEGIN_GROUP(LUACompiler_createVM, "LUA2C");
	UG_DLOG(DID_LUACOMPILER, 2, "parsing " << functionName << "... ");
	try{
		m_f=NULL;
		LUAParserClass parser;
		int ret = parser.parse_luaFunction(functionName);
		if(ret == LUAParserClass::LUAParserError)
		{
			UG_LOG("parsing " << functionName << " failed: reduced LUA parser failed.\n");
			return false;
		}
		if(ret == LUAParserClass::LUAParserIgnore)
		{
			UG_DLOG(DID_LUACOMPILER, 3, "parsing " << functionName << " : Found --LUACompiler:ignore.\n");
			return false;
		}
		//parser.reduce();
	
		string p = PathProvider::get_path(ROOT_PATH) + "/bin/LUACompiler_tmp/";

		if(!DirectoryExists(p))
			CreateDirectory(p);

		fstream out((p + "LUACompiler_output.c").c_str(), fstream::out);
	
		out << "#include <math.h>\n";
		out << "#define true 1\n";
		out << "#define false 0\n";
		parser.createC(out);

		m_iIn = parser.num_in();
		m_iOut = parser.num_out();
		out.close();

		string c1s=string("gcc -fpic -O3 -c ") + p + "LUACompiler_output.c -o " + p + "LUACompiler_output.o";
		UG_DLOG(DID_LUACOMPILER, 2, "compiling line: " << c1s << "\n");
		if(system(c1s.c_str()) != 0)
		{
//			IF_DEBUG(DID_LUACOMPILER, 1)
			if(GetLogAssistant().is_output_process())
			{
				UG_LOG("Error when compiling " << functionName << "\n");
				UG_LOG("compiling line: " << c1s << "\n");
				UG_LOG("--[LUACompiler]-----------------------------------------------\n");
				UG_LOG("created C function from LUA function " << functionName << ":\n");
				UG_LOG(GetFileLines((p+"LUACompiler_output.c").c_str(), 1, -1, true) << "\n");
				UG_LOG("--[LUACompiler]-----------------------------------------------\n");
			}
			return false;
		}
	
		bool bTmpFileSuccess=false;
		m_name = functionName;

		m_pDyn = MakeTmpFile(p+string(functionName), ".dylib", bTmpFileSuccess);

		string c2s=string("gcc -dynamiclib ") + p+"LUACompiler_output.o -o " + m_pDyn.c_str();
		UG_DLOG(DID_LUACOMPILER, 2, "linking line: " << c2s << "\n");
		if(system(c2s.c_str()) != 0)
		{
//			IF_DEBUG(DID_LUACOMPILER, 1)
			if(GetLogAssistant().is_output_process())
			{
				UG_LOG("Error when linking " << functionName << "\n");
				UG_LOG("linking line: " << c2s << "\n");
				UG_LOG("--[LUACompiler]-----------------------------------------------\n");
				UG_LOG("created C function from LUA function " << functionName << ":\n");
				UG_LOG(GetFileLines((p+"LUACompiler_output.c").c_str(), 1, -1, true) << "\n");
				UG_LOG("--[LUACompiler]-----------------------------------------------\n");
			}
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

		if(m_f !=NULL) { UG_DLOG(DID_LUACOMPILER, 2, "OK\n"); }
		else { UG_DLOG(DID_LUACOMPILER, 2, "FAILED\n"); }
		if(m_f !=NULL)
			bInitialized = true;
		return m_f != NULL;
	}
	catch(...)
	{
		UG_DLOG(DID_LUACOMPILER, 1, "exception thrown in LUACompiler::create(" << functionName << ")\n");
		return false;
	}
#else
	UG_THROW("LUA2C not enabled (use LUA2VM).")
#endif

}

bool LUACompiler::createVM(const char *functionName)
{
	PROFILE_BEGIN_GROUP(LUACompiler_createVM, "LUA2VM");
	m_name = functionName;

	LUAParserClass parser;
	try{
		int ret = parser.parse_luaFunction(functionName);
		if(ret == LUAParserClass::LUAParserError)
		{
			UG_LOG("parsing " << functionName << " failed: reduced LUA parser failed.\n");
			return false;
		}
		if(ret == LUAParserClass::LUAParserIgnore)
		{
			UG_DLOG(DID_LUACOMPILER, 3, "parsing " << functionName << " : Found --LUACompiler:ignore.\n");
			return false;
		}

		if(parser.createVM(vm) == false)
		{
			UG_LOG("parsing " << functionName << " failed: create VM failed.\n");
			return false;
		}
		UG_DLOG(DID_LUACOMPILER, 2, "parsing " << functionName << " OK.\n");
	}
	catch(UGError e)
	{
		UG_DLOG(DID_LUACOMPILER, 1, "parsing " << functionName << "... ");
		UG_DLOG(DID_LUACOMPILER, 1, "failed:\n" << e.get_stacktrace() << "\n");
		return false;
	}
	catch(...)
	{
		UG_DLOG(DID_LUACOMPILER, 1, "parsing " << functionName << "... ");
		UG_DLOG(DID_LUACOMPILER, 1, "failed: Exception.\n");
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


LUACompiler::~LUACompiler()
{
    UG_DLOG(DID_LUACOMPILER, 2, "removing " << m_name << "\n");
	if(m_libHandle)
	   	CloseLibrary(m_libHandle);
        
    if(m_pDyn.size() > 0)
    {
		string s = string("rm ") + m_pDyn;
		UG_DLOG(DID_LUACOMPILER, 2, s << "\n");
		system(s.c_str());
	}
}


}
}
