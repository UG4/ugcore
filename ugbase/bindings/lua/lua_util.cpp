/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

//extern libraries
#include <cassert>
#include <cstring>
#include <string>
#include <stack>

// ug libraries
#include "ug.h"

#include "bridge/bridge.h"

#include "common/util/path_provider.h"
#include "common/util/file_util.h"

#include "registry/class_name_provider.h"
#include "registry/class_helper.h"
#include "registry/registry.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "pcl/pcl_util.h"
#endif

#ifdef USE_LUAJIT
#include <lua.hpp>
#else
#include "externals/lua/src/lua.hpp"
#endif


#ifdef UG_DISC
#include "lua_user_data.h"
#endif

#include "bindings_lua.h"
#include "info_commands.h"
#include "lua_stack.h"
#include "lua_debug.h"
#include "lua_util.h"

using namespace std;

namespace ug
{
#ifdef USE_LUA2C
	namespace bridge{
		bool RegisterConverter(Registry &reg, const char* parentGroup);
	}
#endif

namespace script
{

#define PL_COULDNT_FIND (-2)
#define PL_COULDNT_READ (-1)


/**
 *
 * @param filename
 * @param returnedFilename
 * @return true if a file can be found at filename
 */
bool GetNormalFilename(const string &filename, string &returnedFilename)
{
	if(FileExists(filename.c_str())==false) {return false;}
	returnedFilename = filename;
	return true;
}


/**
 * assume UG4_ROOT is the root path of ug4 (containing ugbase, apps, scripts and so on)
 * This function will search a file (in this order)
 * 1.) relative to the current script path.
 * 2.) as a normal, i.e. absolute or relative to the working directory, filename
 * 3.) relative to SCRIPT_PATH (normally UG4_ROOT/scripts)
 * 4.) relative to APPS_PATH (normally UG4_ROOT/apps)
 * 5.) relative to UG4_ROOT
 * @param filename filename to search for
 * @param absoluteFilename absolut filename where the file was found
 * @return true if a file could be found at one of the locations
 */
bool GetAbsoluteUGScriptFilename(const string &filename, string &absoluteFilename)
{
	PROFILE_FUNC();
	return PathProvider::get_filename_relative_to_current_path(filename, absoluteFilename)
			|| GetNormalFilename(filename, absoluteFilename)
			|| PathProvider::get_filename_relative_to_path(SCRIPT_PATH, filename, absoluteFilename)
			|| PathProvider::get_filename_relative_to_path(APPS_PATH, filename, absoluteFilename)
			|| PathProvider::get_filename_relative_to_path(ROOT_PATH, filename, absoluteFilename);
}


std::string GetAbsoluteUGScriptFilenamePaths()
{
	std::stringstream ss;
	ss << "Current Path = " << PathProvider::get_current_path()
			<< ", SCRIPT_PATH = " << PathProvider::get_path(SCRIPT_PATH)
			<< ", APPS_PATH = " << PathProvider::get_path(APPS_PATH)
			<< ", ROOT_PATH = " << PathProvider::get_path(ROOT_PATH);
	return ss.str();
}

/**
 *
 * @param filename			the 'relative' script name. \sa GetAbsoluteUGScriptFilename
 * @param bDistributed	true if loading should be done in parallel
 * @param bThrowOnError		true if errors should be thrown as UGError, else as return false
 * @return true if script could be found and read, false (or UGError) if not
 */
bool LoadUGScript(const char *filename, bool bDistributed, bool bThrowOnError = false)
{
	PROFILE_FUNC();
	string _filename = filename;
	string absoluteFilename = _filename;

	//UG_LOG("LoadUGScript("<<filename<<", "<<(bDistributedLoad?"true":"false") << "\n");
	bool bSuccess = true;
	std::vector<char> file;

#ifdef UG_PARALLEL
	if(pcl::NumProcs() == 1) {bDistributed = false;}
	if(pcl::ProcRank() == 0 || bDistributed==false) {
		bSuccess = GetAbsoluteUGScriptFilename(_filename, absoluteFilename);
	}
	bSuccess = pcl::AllProcsTrue(bSuccess);
#else
	bSuccess = GetAbsoluteUGScriptFilename(_filename, absoluteFilename);
#endif
	if(bSuccess == false)
	{
		if(bThrowOnError)
		{
			UG_THROW("Couldn't find script " << filename << endl
					<< "Search paths: " << GetAbsoluteUGScriptFilenamePaths() << endl);
		}
		return false;
	}

#ifdef UG_PARALLEL
	if(pcl::ParallelReadFile(absoluteFilename, file, true, bDistributed) == false)
#else
	if(ReadFile(absoluteFilename.c_str(), file, true) == false)
#endif
	{
		if(bThrowOnError) { UG_THROW("Couldn't read script " << absoluteFilename << endl); }
		return false;
	}


	// the current script working path is the path of the script
	std::string curPath = PathFromFilename(absoluteFilename);
	PathProvider::push_current_path(curPath);
	{
		// parse and execute the script.
		bSuccess = ParseAndExecuteBuffer(&file[0], (string("@")+absoluteFilename).c_str());
	}
	PathProvider::pop_current_path();
	return bSuccess;
}




bool LoadUGScript_Parallel(const char* filename)
{
	return LoadUGScript(filename, true, true);
}

bool LoadUGScript_Single(const char* filename)
{
	return LoadUGScript(filename, false, true);
}

static bridge::Registry* g_pRegistry = nullptr;

static void UpdateScriptAfterRegistryChange(bridge::Registry* pReg)
{
	PROFILE_FUNC();
	UG_ASSERT(pReg == g_pRegistry, "static g_pRegistry does not match parameter pReg, someone messed up the registries!");
	
//	this can be called since CreateBindings automatically avoids
//	double registration
	bridge::lua::CreateBindings_LUA(GetDefaultLuaState(),
										*pReg);
}

void RegisterDefaultLuaBridge(bridge::Registry* reg, const std::string& grp)
{
	if(reg->functionname_registered("ug_load_script"))
		return;
	
	reg->add_function("ug_load_script", &LoadUGScript_Parallel, "/ug4/lua",
				"success", "", "ONLY IF ALL CORES INVOLVED! Loads and parses a script and returns whether it succeeded.");
	reg->add_function("ug_load_script_single",
			&LoadUGScript_Single, "/ug4/lua",
				"success", "", "Loads and parses a script and returns whether it succeeded.");

	RegisterLuaDebug(*reg);

//	this define makes sure that no methods are referenced that
//	use the algebra, even if no algebra is included.
#ifdef UG_ALGEBRA
//	Register info commands
	RegisterInfoCommands(*reg, grp.c_str());
#ifdef USE_LUA2C
	RegisterConverter(*reg, grp.c_str());
#endif

//	Register user functions
	RegisterLuaUserData(*reg, grp);
#endif
}

static lua_State* theLuaState = nullptr;
lua_State* GetDefaultLuaState()
{
//	if the state has not already been opened then do it now.
	if(!theLuaState)
	{
		PROFILE_BEGIN(CreateLUARegistry);
		if(!g_pRegistry){
		//	store a pointer to the registry and avoid multiple callback registration
			g_pRegistry = &bridge::GetUGRegistry();
			g_pRegistry->add_callback(UpdateScriptAfterRegistryChange);
		}
		
	//	open a lua state
		theLuaState = luaL_newstate();
#ifdef USE_LUAJIT
		UG_ASSERT(theLuaState!=nullptr, "FATAL ERROR: Not enough memory for lua?")
#endif

	//	open standard libs
		luaL_openlibs(theLuaState);

	//	make metatables available in lua script
		lua_register(theLuaState, "ug_get_metatable", UGGetMetatable );

	//	make base class check available in lua script
		lua_register(theLuaState, "ug_is_base_class", UGIsBaseClass);

	//	make dim check available in lua script
		lua_register(theLuaState, "ug_dim_compiled", UGDimCompiled);

	//	make algebra check available in lua script
		lua_register(theLuaState, "ug_algebra_compiled", UGAlgebraCompiled);

	//	make class name available in lua script
		lua_register(theLuaState, "ug_class_name", UGGetClassName);
		lua_register(theLuaState, "ug_class_group", UGGetClassGroup);

	//	create lua bindings for registered functions and objects
		bridge::lua::CreateBindings_LUA(theLuaState, *g_pRegistry);
	}
	
	return theLuaState;
}

/// calls lua_close, which calls delete for all lua objects
void ReleaseDefaultLuaState()
{
	if(theLuaState != nullptr)
	{
		lua_close(theLuaState);
		theLuaState = nullptr;
	}
	FinalizeLUADebug();
}


/// error function to be used for lua_pcall
int luaCallStackError( lua_State *L )
{
	UG_LOG("LUA-ERROR! Call stack:\n");
    bridge::LuaStackTrace(0);
    return 1;
}


bool ParseAndExecuteBuffer(const char* buffer, const char *bufferName)
{
	PROFILE_FUNC();
	lua_State* L = GetDefaultLuaState();

	lua_pushcfunction(L, luaCallStackError);

	PROFILE_BEGIN(luaL_loadbuffer);
	int error = luaL_loadbuffer(L, buffer, strlen(buffer), bufferName);
	PROFILE_END_(luaL_loadbuffer);

	if(error == 0)
	{
		PROFILE_BEGIN(lua_pcall);
		error = lua_pcall(L, 0, 0, -2);
	}

	if(error)
	{
		const string msg = lua_tostring(L, -1);
		lua_pop(L, 1);
		if(msg.find("__UG__LUA__EMPTY__MSG__") == string::npos)
			throw LuaError(msg.c_str());
		throw LuaError();
	}
	return true;
}


static void GetToStringFromStack(lua_State *L, std::stringstream &ss)
{
	const int nArgs = lua_gettop(L);
	lua_getglobal(L, "tostring");

	for(int i = 1; i<=nArgs; i++)
	{
		lua_pushvalue(L, -1);
		lua_pushvalue(L, i);
		lua_call(L, 1, 1);
		const char *s = lua_tostring(L, -1);
		if(s) ss << s;
		lua_pop(L, 1);
	}
	lua_pop(L,1);
}

/// UGLuaPrint. Redirects LUA prints to UG_LOG
int UGLuaPrint(lua_State *L)
{
	std::stringstream ss;
	GetToStringFromStack(L, ss);
	ss << "\n";
	UG_LOG(ss.str());
	return 0;
}

/// UGLuaPrint. Redirects LUA prints to UG_LOG
int UGLuaPrintAllProcs(lua_State *L)
{
	std::stringstream ss;
	GetToStringFromStack(L, ss);
	ss << "\n";
	UG_LOG_ALL_PROCS(ss.str());
	return 0;
}


/// UGLuaWrite. Redirects LUA prints to UG_LOG without adding newline at the end
int UGLuaWrite(lua_State *L)
{
	std::stringstream ss;
	GetToStringFromStack(L, ss);
	UG_LOG(ss.str());
	GetLogAssistant().flush();
	return 0;
}

/// UGLuaErrLog.
int UGLuaErrLog(lua_State *L)
{
	std::stringstream ss;
	GetToStringFromStack(L, ss);
	ss << "\n";
	UG_ERR_LOG(ss.str());
	GetLogAssistant().flush();
	return 0;
}

void RegisterStdLUAFunctions(lua_State *L)
{
	lua_register(L, "print", UGLuaPrint );
	lua_register(L, "print_all", UGLuaPrintAllProcs );
	lua_register(L, "write", UGLuaWrite );
	lua_register(L, "err_log", UGLuaErrLog);
}

int UGGetMetatable(lua_State *L)
{
	const char* name = lua_tostring(L, -1);
	if(name){
		luaL_getmetatable(L, name);
		return 1;
	}
	lua_pushnil(L);
	return 1;
}

int UGGetClassName(lua_State *L)
{
	if(lua_getmetatable(L, -1) != 0)
	{
		lua_pushstring(L, "class_name_node");
		lua_rawget(L, -2);
		const auto* classNameNode = static_cast<const bridge::ClassNameNode *>(lua_touserdata(L, -1));
		lua_pop(L, 2);

		if(classNameNode)
		{
			lua_pushstring(L, classNameNode->name().c_str());
			return 1;
		}
		else UG_THROW("In UGGetClassName: Something wrong with ClassNameNode.");
	}

	lua_pushstring(L, "");
	return 1;
}

int UGGetClassGroup(lua_State *L)
{
	if(lua_getmetatable(L, -1) != 0)
	{
		lua_pushstring(L, "class_name_node");
		lua_rawget(L, -2);
		const auto* classNameNode = static_cast<const bridge::ClassNameNode *>(lua_touserdata(L, -1));
		lua_pop(L, 2);

		if(classNameNode)
		{
			const bridge::Registry &reg = bridge::GetUGRegistry();
			const bridge::IExportedClass *c = reg.get_class(classNameNode->name());
			lua_pushstring(L, c->group().c_str());
			return 1;
		}
		else UG_THROW("In UGGetClassGroup: Something wrong with ClassNameNode.");
	}

	lua_pushstring(L, "");
	return 1;
}


int UGIsBaseClass(lua_State *L)
{
	const char* derivClass = lua_tostring(L, -1);
	const char* baseClass = lua_tostring(L, -2);

	if(!derivClass || !baseClass)
		UG_THROW("In UGIsBaseClass: invalid class names passed as arguments.");

	luaL_getmetatable(L, derivClass);

	if(!lua_isnil(L, -1))
	{
		lua_pushstring(L, "class_name_node");
		lua_rawget(L, -2);
		const auto* classNameNode = static_cast<const bridge::ClassNameNode *>(lua_touserdata(L, -1));
		lua_pop(L, 2);

		if(classNameNode)
		{
			if(!classNameNode->empty())
			{
				if(ClassNameTreeContains(*classNameNode, baseClass))
				{
					lua_pushboolean(L, true);
					return 1;
				}
			}
			lua_pushboolean(L, false);
			return 1;
		}
		else UG_THROW("In UGIsBaseClasse: Something wrong with ClassNameNode.");
	}
	else{
		UG_THROW("In UGIsBaseClass: Cannot find metatable for: "<< derivClass);
	}
	return 0;
}

void SetLuaUGArgs(lua_State* L, int argc, char* argv[])
{
	lua_newtable(L);
	int ugargc=0;
	for(int i = 0; i < argc; ++i){
		lua_pushnumber(L, ++ugargc); //	push the index to the table
		lua_pushstring(L, argv[i]); //	push the value to the table
		lua_settable(L, -3); //	create the entry
	}

	lua_setglobal(L, "ugargv");//	set the tables name
	lua_pushnumber(L, ugargc);
	lua_setglobal(L, "ugargc");
}


int UGDimCompiled(lua_State *L)
{
	int dim = lua_tointeger(L, -1);
#if UG_DIM_1
	if(dim == 1){ lua_pushboolean(L, true); return 1;}
#endif
#if UG_DIM_2
	if(dim == 2){ lua_pushboolean(L, true); return 1;}
#endif
#if UG_DIM_3
	if(dim == 3){ lua_pushboolean(L, true); return 1;}
#endif
	lua_pushboolean(L, false); return 1;
}

int UGAlgebraCompiled(lua_State *L)
{
	std::string name = lua_tostring(L, -1);
#if UG_CPU_1
	if(name == "CPU1"){ lua_pushboolean(L, true); return 1;}
#endif
#if UG_CPU_2
	if(name == "CPU2"){ lua_pushboolean(L, true); return 1;}
#endif
#if UG_CPU_3
	if(name == "CPU3"){ lua_pushboolean(L, true); return 1;}
#endif
#if UG_CPU_4
	if(name == "CPU4"){ lua_pushboolean(L, true); return 1;}
#endif
#if UG_CPU_5
	if(name == "CPU5"){ lua_pushboolean(L, true); return 1;}
#endif
#if UG_CPU_6
	if(name == "CPU6"){ lua_pushboolean(L, true); return 1;}
#endif
#if UG_CPU_VAR
	if(name == "CPUVAR"){ lua_pushboolean(L, true); return 1;}
#endif
	lua_pushboolean(L, false); return 1;
}
}
}