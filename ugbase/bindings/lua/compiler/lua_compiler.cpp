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
#include "common/util/file_util.h"
#include "common/util/path_provider.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "lua_compiler.h"
#include "lua_compiler_debug.h"
#include "common/profiler/profiler.h"
#include "vm.h"
using namespace std;

namespace ug{

extern bool useLua2VM;
extern bool useLuaCompiler;

DebugID DID_LUACOMPILER("LUACompiler");

namespace bridge {
    

bool LUACompiler::create(const char *functionName, LuaFunctionHandle* pHandle)
{
	if(useLua2VM)
		return createVM(functionName, pHandle);
	else
		return createC(functionName, pHandle);
}

bool LUACompiler::createC(const char *functionName, LuaFunctionHandle* pHandle)
{
#ifdef USE_LUA2C
	PROFILE_BEGIN_GROUP(LUACompiler_createVM, "LUA2C");
	UG_DLOG(DID_LUACOMPILER, 1, "LUA2C: parsing " << functionName << "... ");
	try{
		m_f=nullptr;
		LUAParserClass parser;
		int ret = 0;
		if(pHandle == nullptr){
			ret = parser.parse_luaFunction(functionName);
		} else {
			ret = parser.parse_luaFunction(*pHandle);
		}
		if(ret == LUAParserClass::LUAParserError)
		{
			UG_DLOG(DID_LUACOMPILER, 1, "failed: reduced LUA parser failed.\n");
			return false;
		}
		if(ret == LUAParserClass::LUAParserIgnore)
		{
			UG_DLOG(DID_LUACOMPILER, 1, "Found --LUACompiler:ignore. Ignoring this function.\n");
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

		ret = parser.createC(out);
		if(ret != LUAParserClass::LUAParserOK)
		{
			UG_DLOG(DID_LUACOMPILER, 1, "some problem when generating C code.\n");
			return false;
		}

		m_iIn = parser.num_in();
		m_iOut = parser.num_out();
		out.close();

		UG_DLOG(DID_LUACOMPILER, 5, GetFileLines((p+"LUACompiler_output.c").c_str(), 1, -1, true) << "\n");


		string c1s=string("gcc -fpic -O3 -c ") + p + "LUACompiler_output.c -o " + p + "LUACompiler_output.o";
		UG_DLOG(DID_LUACOMPILER, 2, "compiling line: " << c1s << "\n");
		if(system(c1s.c_str()) != 0)
		{
//			IF_DEBUG(DID_LUACOMPILER, 1)
			if(GetLogAssistant().is_output_process())
			{
				UG_LOG("\nLUA2C: Error when compiling " << functionName << "\n");
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

		if(GetLogAssistant().is_output_process())
		{	UG_DLOG(DID_LUACOMPILER, 2, "linking line: " << c2s << "\n"); }

		if(system(c2s.c_str()) != 0)
		{
//			IF_DEBUG(DID_LUACOMPILER, 1)
			if(GetLogAssistant().is_output_process())
			{
				UG_LOG("\nLUA2C: Error when linking " << functionName << "\n");
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
			UG_LOG("\nLUA2C: Error when opening library for function " << functionName << "\n");
			UG_LOG("Error is " << error << "\n");
			return false;
		}
		m_f = (LUA2C_Function) GetLibraryProcedure(m_libHandle, functionName);

		if(m_f !=nullptr) { UG_DLOG(DID_LUACOMPILER, 1, "OK\n"); }
		else { UG_DLOG(DID_LUACOMPILER, 1, "FAILED\n"); }
		if(m_f !=nullptr)
			bInitialized = true;
		return m_f != nullptr;
	}
	catch(...)
	{
		UG_DLOG(DID_LUACOMPILER, 1, "LUA2C: exception thrown in LUACompiler::create(" << functionName << ")\n");
		return false;
	}
#else
	UG_THROW("LUA2C not enabled (use LUA2VM).")
#endif

}

bool LUACompiler::createVM(const char *functionName, LuaFunctionHandle* pHandle)
{
	PROFILE_BEGIN_GROUP(LUACompiler_createVM, "LUA2VM");
	m_name = functionName;


	LUAParserClass parser;
	try
	{
		int ret = 0;
		if(pHandle == nullptr){
			parser.parse_luaFunction(functionName);
		} else {
			parser.parse_luaFunction(*pHandle);
		}
		if(vm != nullptr) delete vm;
		vm = new VMAdd;
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

		if(parser.createVM(*vm) == false)
		{
			UG_LOG("parsing " << functionName << " failed: create VM failed.\n");
			return false;
		}

		IF_DEBUG(DID_LUACOMPILER, 5)
		{	vm->print(); }

		UG_DLOG(DID_LUACOMPILER, 1, "LUA2VM: parsing " << functionName << " OK.\n");
	}
	catch(UGError e)
	{
		UG_DLOG(DID_LUACOMPILER, 1, "LUA2VM: parsing " << functionName << "... ");
		UG_DLOG(DID_LUACOMPILER, 1, "failed:\n" << e.get_stacktrace() << "\n");
		return false;
	}
	catch(...)
	{
		UG_DLOG(DID_LUACOMPILER, 1, "LUA2VM: parsing " << functionName << "... ");
		UG_DLOG(DID_LUACOMPILER, 1, "failed: Exception.\n");
		return false;
	}
	//UG_LOG(" ok.\n");
	m_iIn = vm->num_in();
	m_iOut = vm->num_out();
	bInitialized = true;
	bVM = true;
	return true;
}


LUACompiler::~LUACompiler()
{
	if(vm != nullptr) delete vm;

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

bool LUACompiler::call(double *ret, const double *in) const
{
	if(bVM)
	{
		const_cast<LUACompiler*>(this)->vm->execute(ret, in);
		return true;
	}
	else
	{
		UG_ASSERT(m_f != nullptr, "function " << m_name << " not valid");
		m_f(ret, in);
		return true;
	}
}


}
}
