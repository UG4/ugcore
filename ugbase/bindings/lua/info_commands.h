/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include "externals/lua/lua.h"
#include "registry/registry.h"
#include <string>

#ifndef __H__UG_SCRIPT__INFO_COMMANDS__
#define __H__UG_SCRIPT__INFO_COMMANDS__

#include "common/ug_config.h"

namespace ug
{
namespace bridge
{

///	registers info commands TypeInfo, ClassUsage and others
UG_API bool RegisterInfoCommands(bridge::Registry &reg, const char* grp);

/**
 * \param L 		the lua state
 * \param index		the index of the lua object
 * \return the corresponing ClassNameNode
 * \sa ClassNameNode
 */
UG_API const ClassNameNode* GetClassNameNode(lua_State *L, int index);

/**
 * \param L 		the lua state
 * \param index		the index of the lua object
 * \return a std::vector of class names
 */
UG_API const std::vector<const char*> *GetClassNames(lua_State *L, int index);

/**
 * \param L 		the lua state
 * \param name		the name of the lua object
 * \return the class names of the object or nullptr if not found
 */
UG_API const std::vector<const char*> *GetClassNames(lua_State* L, const char *name);


/**
 * \brief prints out information for a method of a class
 * \param L 			the lua state
 * \param index			index of the class object on the lua stack
 * \param thefunc		method to print
 * \return 0
 */
UG_API std::string LuaClassMethodInfo(lua_State *L, int index, const ExportedMethod &thefunc);



/**
 * \brief prints out information for a lua function (a function defined in lua script)
 * \param L 			the lua state
 * \param bComplete		if complete, we print the source code of the function
 * \return 0
 */
UG_API std::string FunctionInfo(lua_State *L, bool bComplete, const char *functionName=nullptr);

UG_API int UGTypeInfo(const char *p);

UG_API bool GetLuaNamespace(lua_State* L, std::string name);

/**
 * \brief returns an integer to a lua-variable.
 * \param L				the lua state
 * \param name			name of the variable. namespaces are possible (like math.pi)
 * \param notAvailable	return value if variable was not found
  */
UG_API double LuaGetNumber(lua_State *L, const char *name, double notAvailable);

/**
 * \brief returns the string to a lua-variable.
 * \param L				the lua state
 * \param name			name of the variable. namespaces are possible (like math.pi)
 * \param notAvailable	return value if variable was not found
  */
UG_API std::string LuaGetString(lua_State *L, const char *name, const char *notAvailable);

/**
 * \brief returns a boolean to a lua-variable.
 * \param L				the lua state
 * \param name			name of the variable. namespaces are possible (like math.pi)
 * \param notAvailable	return value if variable was not found
  */
UG_API bool LuaGetBoolean(lua_State *L, const char *name, bool notAvailable);

/**
 * prints the current lua line
 * @param L
 */
UG_API std::string LuaCurrentLine(lua_State* L);

/**
 * searches the lua stack for the first valid line info
 * @param L			the lua state
 * @param entry 	returned lua_Debug structure
 */
UG_API void LuaGetLastLine(lua_State* L, lua_Debug entry);

/**
 * \brief prints the source of a lua script function which is on top of the stack
 * \param L		the lua state
 */
UG_API std::string LuaGetScriptFunctionString(lua_State *L, int index=-1);

/**
 * \brief prints the source of a lua script function which is on top of the stack
 * \param L			the lua state
 * \param index 	lua stack index of the table
 */
UG_API void LuaPrintTable(lua_State *L, size_t iSpace, int index=-1);
UG_API bool ClassInstantiations(const char *classname);


/**
 * returns a String describing the content of the lua stack at a given index
 * and all types it is compatible with (for ex. "2" is string and number)
 * @param L			the lua state
 * @param index		lua stack index to get type string of
 * @return			type string
 */
UG_API std::string GetLuaTypeString(lua_State* L, int index);

///
/**
 * @brief  prints information about lua's call stack (file:line source).
 * @param L the lua state
 * @param backtraceLevel number of calls to display (<= 0 -> all)
 */
UG_API std::string LuaStackTraceString(lua_State* L, int fromLevel=0, int toLevel=-1);

UG_API std::string LuaStackTraceString();

/// prints information about lua's call stack (file:line source).
UG_API void LuaStackTrace(int fromlevel);


/// returns the current file and line-number and line-content ( \sa LuaStackTrace ).
UG_API std::string GetLuaFileAndLine(lua_State* L);

/// returns the current file and line-number
UG_API std::string GetLuaFileAndLineNumber(lua_State* L);

/// returns the current line-content
UG_API std::string GetLuaLine(lua_State* L);

/// returns file and line of defined script function
std::string GetLUAScriptFunctionDefined(const char *functionName);

}
}
						  
#endif // __H__UG_SCRIPT__INFO_COMMANDS__
