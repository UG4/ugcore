/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef UG_BASE_BINDINGS_LUA_LUA_UTIL_H
#define UG_BASE_BINDINGS_LUA_LUA_UTIL_H

#include <string>

#ifdef USE_LUAJIT
#include <lua.hpp>
#else
// default lua
extern "C" {
#include "externals/lua/src/lua.hpp" // for apps/unit_test
}
#endif

#include "common/common.h"
#include "registry/registry.h"



namespace ug {
namespace script {

///	Error class thrown if an error occurs during parsing.
class LuaError final : public UGError
{
	public:
		explicit LuaError(const char* msg) : UGError(msg), bShowMsg(true)	{}
		LuaError() : UGError(""), bShowMsg(false)	{}

		bool show_msg() const {return bShowMsg;}

	protected:
		bool bShowMsg;
};


///	loads and parses a file. Several paths are tried if the file is not found.
/**	Throws an instance of LuaError, if a parse error occurs.
 * This method first tries to load the file specified with filename relative
 * to the path of the currently parsed file (if LoadUGScript is called from
 * within a load-script). If this failed, the file is tried to be loaded
 * with the raw specified filename. If this fails too, the method tries to
 * load the file from ugs scripting directory.
 *
 * Note that this method pushes the path of the currently parsed script to
 * PathProvider when parsing starts, and pops it when parsing is done.
 * \param filename The filename for the script to be loaded. may be relative to ug4/apps/
 * \param bDistributed if true, loads the script on core 0 and distributes
 * \param bThrowOnError if true throws on error, else returns false
 * the script to all other cores via pcl::broadcast . Use this whenever
 * possible when doing parallel work otherwise this routine can take a long time
 * (ignored if UG_PARALLEL not defined)
 */
UG_API bool LoadUGScript(const char* filename, bool bDistributed, bool bThrowOnError);
/// calls LoadUGScript with bDistributed=true
UG_API bool LoadUGScript_Parallel(const char* filename);
/// calls LoadUGScript with bDistributed=false . Avoid. (\sa LoadUGScript)
UG_API bool LoadUGScript_Single(const char* filename);

/// registers lua only functionality at the registry
UG_API void RegisterDefaultLuaBridge(bridge::Registry* reg, const std::string& grp = "/ug4");

///	returns the default lua state
/**	When called for the first time, or after ReleaseDefaultLuaState,
 * a new state is created and the methods and classes in ugs default registry
 * (ug::bridge::GetUGRegistry) are registered. Furthermore a callback
 * is registered, which registers new methods whenever
 * Registry::registry_changed() is called on the default registry.*/
UG_API lua_State* GetDefaultLuaState();

/// Releases the lua-state returned by GetDefaultLuaState().
/**	This method is useful, if you want to restart scripting from scratch.*/
UG_API void ReleaseDefaultLuaState();

/**
 * Parses the content of buffer and executes it in the default lua state
 * @param buffer		the buffer to be executed
 * @param bufferName	name of the buffer (for error messages)
 * @return				true on success, otherwise throw(LuaError)
 * Throws an instance of LuaError, if a parse error occurs.
 */
UG_API bool ParseAndExecuteBuffer(const char* buffer, const char *bufferName="buffer");

/// UGLuaPrint. Redirects LUA prints to UG_LOG
UG_API int UGLuaPrint(lua_State *L);

/// UGLuaPrintAllProcs. Redirects LUA prints to UG_LOG_ALL_PROCS
int UGLuaPrintAllProcs(lua_State *L);

/// UGLuaWrite. prints LUA output to UG_LOG without adding std::endl automatically
UG_API int UGLuaWrite(lua_State *L);

///	Returns the metatable for the given class
UG_API int UGGetMetatable(lua_State *L);

///	Returns if a class contains a base class
UG_API int UGIsBaseClass(lua_State *L);

/// Returns if dimension is compiled into binary
UG_API int UGDimCompiled(lua_State *L);

/// Returns if dimension is compiled into binary
UG_API int UGAlgebraCompiled(lua_State *L);

/// Returns type of a userdata as string
UG_API int UGGetClassName(lua_State *L);

/// Returns classgroup of a userdata as string
UG_API int UGGetClassGroup(lua_State *L);

/**
 * create ugargc and ugargv in lua
 * we want to forward argc and argv to the lua-environment.
 * we'll create a table for that.
 * @param L					the LUA state
 * @param argc				the commandline argc
 * @param argv				the commandline argv
 */
UG_API void SetLuaUGArgs(lua_State* L, int argc, char* argv[]);

/// register functions like print and write directly to LUA (not using the ug registry)
UG_API void RegisterStdLUAFunctions(lua_State *L);

/**
 * searches for the filename
 * - relative to current script
 * - as absolute filename
 * - in PathProvider::get_path(SCRIPT_PATH) (ug4/scripts)
 * - in PathProvider::get_path(APPS_PATH) (ug4/apps)
 * - in PathProvider::get_path(ROOT_PATH) (ug4)
 * \param filename in: relative filename to paths above. out: absolute filename (if found)
 * \param absoluteFilename returns the absolut filename where the file was found
 * \return true if found, else false
 */
UG_API bool GetAbsoluteUGScriptFilename(const std::string &filename, std::string &absoluteFilename);

}
}

#endif