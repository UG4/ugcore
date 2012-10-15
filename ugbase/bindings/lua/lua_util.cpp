// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d20

//extern libraries
#include <cassert>
#include <cstring>
#include <string>
#include <stack>


// ug libraries
#include "ug.h"
#include "lua_util.h"
#include "common/util/path_provider.h"
#include "common/util/file_util.h"
#include "bindings_lua.h"
#include "bridge/bridge.h"
#include "registry/class_helper.h"
#include "info_commands.h"
#include "lua_user_data.h"
#include "registry/class_name_provider.h"
#include "registry/registry.h"
#include "lua_debug.h"

#include "common/util/binary_buffer.h"


#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "pcl/pcl_util.h"
#endif


using namespace std;

namespace ug
{

namespace script
{


#define PL_COULDNT_FIND -2
#define PL_COULDNT_READ -1



bool GetAbsoluteFilename(const string &relativeFilename, string &absoluteFilename)
{
	if(FileExists(relativeFilename.c_str())==false) return false;
	absoluteFilename = relativeFilename;
	return true;
}


bool GetAbsoluteUGScriptFilename(const string &filename, string &absoluteFilename)
{
	PROFILE_FUNC();
	return PathProvider::get_filename_relative_to_current_path(filename, absoluteFilename)
			|| GetAbsoluteFilename(filename, absoluteFilename)
			|| PathProvider::get_filename_relative_to_path(SCRIPT_PATH, filename, absoluteFilename)
			|| PathProvider::get_filename_relative_to_path(APPS_PATH, filename, absoluteFilename)
			|| PathProvider::get_filename_relative_to_path(ROOT_PATH, filename, absoluteFilename);
}



bool LoadUGScript(const char *_filename, bool bDistributedLoad)
{
	PROFILE_FUNC();
	string filename=_filename;
	string absoluteFilename=filename;

	//UG_LOG("LoadUGScript("<<filename<<", "<<(bDistributedLoad?"true":"false") << "\n");
	bool bSuccess = true;
	std::vector<char> file;

#ifdef UG_PARALLEL
	if(pcl::GetNumProcesses() == 1) bDistributedLoad = false;
	if(pcl::GetProcRank() == 0 || bDistributedLoad==false)
		bSuccess = GetAbsoluteUGScriptFilename(filename, absoluteFilename);
	bSuccess = pcl::AllProcsTrue(bSuccess);
#else
	bSuccess = GetAbsoluteUGScriptFilename(filename, absoluteFilename);
#endif
	if(bSuccess == false)
	{
		UG_LOG("Couldn't find script " << _filename << endl);
		return false;
	}
#ifdef UG_PARALLEL
	if(pcl::ParallelReadFile(absoluteFilename, file, true, bDistributedLoad) == false)
#else
	if(ReadFile(absoluteFilename.c_str(), file, true) == false)
#endif
	{
		UG_LOG("Couldn't read script " << absoluteFilename << endl);
		return false;
	}

	std::string curPath = PathFromFilename(absoluteFilename);
	PathProvider::push_current_path(curPath);
	bool success = ParseBuffer(&file[0], (string("@")+absoluteFilename).c_str());
	//bool success = ParseFile(absoluteFilename.c_str());
	PathProvider::pop_current_path();

	return success;
}




bool LoadUGScript_Parallel(const char* filename)
{
	return LoadUGScript(filename, true);
}
bool LoadUGScript_Single(const char* filename)
{
	return LoadUGScript(filename, false);
}

static ug::bridge::Registry* g_pRegistry = NULL;

static void UpdateScriptAfterRegistryChange(ug::bridge::Registry* pReg)
{
	PROFILE_FUNC();
	UG_ASSERT(pReg == g_pRegistry, "static g_pRegistry does not match parameter pReg, someone messed up the registries!");
	
//	this can be called since CreateBindings automatically avoids
//	double registration
	ug::bridge::lua::CreateBindings_LUA(GetDefaultLuaState(),
										*pReg);
}


void RegisterDefaultLuaBridge(ug::bridge::Registry* reg, std::string grp)
{

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

//	Register user functions
	RegisterLuaUserData(*reg, grp);
	#endif
}

static lua_State* theLuaState = NULL;
lua_State* GetDefaultLuaState()
{
//	if the state has not already been opened then do it now.
	if(!theLuaState)
	{
		PROFILE_BEGIN(CreateLUARegistry);
		if(!g_pRegistry){
		//	store a pointer to the registry and avoid multiple callback registration
			g_pRegistry = &ug::bridge::GetUGRegistry();
			g_pRegistry->add_callback(UpdateScriptAfterRegistryChange);
		}
		
	//	open a lua state
		theLuaState = lua_open();
	//	open standard libs
		luaL_openlibs(theLuaState);

	//	make metatables available in lua script
		lua_register(theLuaState, "ug_get_metatable", UGGetMetatable );

	//	make base class check available in lua script
		lua_register(theLuaState, "ug_is_base_class", UGIsBaseClass);

	//	make dim check available in lua script
		lua_register(theLuaState, "ug_dim_compiled", UGDimCompiled);

	//	make class name available in lua script
		lua_register(theLuaState, "ug_class_name", UGGetClassName);

	//	create lua bindings for registered functions and objects
		ug::bridge::lua::CreateBindings_LUA(theLuaState, *g_pRegistry);
		PROFILE_END();
	}
	
	return theLuaState;
}

/// calls lua_close, which calls delete for all lua objects
void ReleaseDefaultLuaState()
{
	if(theLuaState != NULL)
	{
		lua_close(theLuaState);
		theLuaState = NULL;
	}
	FinalizeLUADebug();
	return;
}


/// error function to be used for lua_pcall
int luaCallStackError( lua_State *L )
{
	//UG_LOG("Error: " << lua_tostring(L, -1) << ". ");
    //UG_LOG("call stack:\n"); ug::bridge::LuaStackTrace(L);
    return 1;
}

/**
 * Parses the content of buffer and executes it in the default lua state
 * @param buffer		the buffer to be executed
 * @param bufferName	name of the buffer (for error messages)
 * @return				true on success, otherwise throw(LuaError)
 */
bool ParseBuffer(const char* buffer, const char *bufferName)
{
	PROFILE_BEGIN(ParseBuffer);
	lua_State* L = GetDefaultLuaState();

	lua_pushcfunction(L, luaCallStackError);

	PROFILE_BEGIN(luaL_loadbuffer);
	int error = luaL_loadbuffer(L, buffer, strlen(buffer), bufferName);
	PROFILE_END();

	if(error == 0)
	{
		PROFILE_BEGIN(lua_pcall);
		error = lua_pcall(L, 0, 0, -2);
	}

	if(error)
	{
//		LOG("PARSE-ERROR: " << lua_tostring(L, -1) << endl);
		string msg = lua_tostring(L, -1);
		lua_pop(L, 1);
		//ug::bridge::LuaStackTrace(L);
		throw(LuaError(msg.c_str()));
//		return false;
	}

	PROFILE_END();
	return true;

}

/**
 * Parses the content of the file "filename" and executes it in the default lua state
 * \note don't use this function on all cores (i/o limited).
 * @param filename		the buffer to be executed
 * @return				true on success, otherwise throw(LuaError)
 */
bool ParseFile(const char* filename)
{
	lua_State* L = GetDefaultLuaState();

	lua_pushcfunction(L, luaCallStackError);

	PROFILE_BEGIN(luaL_loadfile);
	int error = luaL_loadfile(L, filename);
	PROFILE_END();

	if(error == 0)
	{
		PROFILE_BEGIN(lua_pcall);
		error = lua_pcall(L, 0, 0, -2);
		PROFILE_END();
	}

	if(error)
	{
		//LOG("PARSE-ERROR in parse_file(" << filename << "): " << lua_tostring(L, -1) << endl);
		string msg = lua_tostring(L, -1);
		lua_pop(L, 1);
		//ug::bridge::LuaStackTrace(L);
		throw(LuaError(msg.c_str()));
	}
	return true;
}

/// UGLuaPrint. Redirects LUA prints to UG_LOG
int UGLuaPrint(lua_State *L)
{
	int nArgs = lua_gettop(L);
	int i;
	lua_getglobal(L, "tostring");

	for(i=1; i<=nArgs; i++)
	{
		lua_pushvalue(L, -1);
		lua_pushvalue(L, i);
		lua_call(L, 1, 1);
		const char *s = lua_tostring(L, -1);
		if(s) UG_LOG(s);
		lua_pop(L, 1);
	}
	UG_LOG(endl);
	lua_pop(L,1);
	return 0;
}

/// UGLuaWrite. Redirects LUA prints to UG_LOG without adding newline at the end
int UGLuaWrite(lua_State *L)
{
	PROFILE_FUNC();
	int nArgs = lua_gettop(L);
	int i;
	lua_getglobal(L, "tostring");

	for(i=1; i<=nArgs; i++)
	{
		lua_pushvalue(L, -1);
		lua_pushvalue(L, i);
		lua_call(L, 1, 1);
		const char *s = lua_tostring(L, -1);
		if(s) UG_LOG(s);
		lua_pop(L, 1);
	}
	lua_pop(L,1);
	return 0;
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
		const ug::bridge::ClassNameNode* classNameNode = (const ug::bridge::ClassNameNode*) lua_touserdata(L, -1);
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
		const ug::bridge::ClassNameNode* classNameNode = (const ug::bridge::ClassNameNode*) lua_touserdata(L, -1);
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
}}//	end of namespace
