
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
UG_API bool RegisterInfoCommands(bridge::Registry &reg, const char* parentGroup = "/ug4");

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
 * \return the class names of the object or NULL if not found
 */
UG_API const std::vector<const char*> *GetClassNames(lua_State* L, const char *name);


/**
 * \brief prints out information for a method of a class
 * \param L 			the lua state
 * \param index			index of the class object on the lua stack
 * \param thefunc		method to print
 * \return 0
 */
UG_API void PrintLuaClassMethodInfo(lua_State *L, int index, const ExportedMethod &thefunc);



/**
 * \brief prints out information for a lua function (a function defined in lua script)
 * \param L 			the lua state
 * \param bComplete		if complete, we print the source code of the function
 * \return 0
 */
UG_API int PrintFunctionInfo(lua_State *L, bool bComplete);

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
UG_API void LuaPrintCurrentLine(lua_State* L);

/**
 * searches the lua stack for the first valid line info
 * @param L			the lua state
 * @param entry 	returned lua_Debug structure
 */
UG_API void LuaGetLastLine(lua_State* L, lua_Debug entry);

/**
 * \brief function to get a line of a file
 * \param filename				file name
 * \param line
 * \return the line of the file.
 */
std::string GetFileLine(const char *filename, size_t line);

/**
 * \brief function to get some lines of a file
 * \param filename				file name
 * \param fromline
 * \param toline
 * \param includeLineNumbers	if true, add the line number in front of each line and a tab.
 *
 * \return lines fromline to toline of file filename.
 */
std::string GetFileLines(const char *filename, size_t fromline, size_t toline, bool includeLineNumbers=false);


/**
 * \brief prints the source of a lua script function which is on top of the stack
 * \param L		the lua state
 */
void PrintLuaScriptFunction(lua_State *L, int index=-1);

/**
 * \brief prints the source of a lua script function which is on top of the stack
 * \param L			the lua state
 * \param index 	lua stack index of the table
 */
void LuaPrintTable(lua_State *L, size_t iSpace, int index=-1);
bool ClassInstantiations(const char *classname);


/**
 * returns a String describing the content of the lua stack at a given index
 * and all types it is compatible with (for ex. "2" is string and number)
 * @param L			the lua state
 * @param index		lua stack index to get type string of
 * @return			type string
 */
std::string GetLuaTypeString(lua_State* L, int index);

/// prints information about lua's call stack (file:line source).
void LuaStackTrace(lua_State* L);

/// returns the current file and line ( \sa LuaStackTrace ).
std::string GetLuaFileAndLine(lua_State* L);

}
}
						  
#endif // __H__UG_SCRIPT__INFO_COMMANDS__
