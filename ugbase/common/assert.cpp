/*
 * assert.cpp
 *
 *  Created on: 05.07.2012
 *      Author: Martin Rupp
 */

#include "common/log.h"
#include "common/profiler/profiler.h"
#include "assert.h"
#include <stdlib.h>

#ifdef UG_POSIX
#include "util/demangle.h"
#include <execinfo.h>
#endif

#include <string>
#include <sstream>

#ifdef UG_FOR_LUA
#include "bindings/lua/info_commands.h"
#include <map>
namespace ug{namespace script{bool IsLUADebug();}}
#endif

#ifdef UG_PARALLEL
#include "pcl/pcl_base.h"
#endif

using namespace std;

void lua_backtrace()
{
#ifdef UG_FOR_LUA
	bool bDoLUAStackTrace=true;
#ifdef UG_POSIX
#ifdef UG_PROFILER_SHINY
	if(ug::script::IsLUADebug()) bDoLUAStackTrace = false;
#endif
#endif
	if(bDoLUAStackTrace)
	{
		UG_LOG("--------------- LUA Stack Trace: --------------\n");
		UG_LOG(ug::bridge::LuaStackTraceString());
		UG_LOG("\n");
	}
#endif
}

void shiny_backtrace()
{
#ifdef UG_POSIX
#ifdef UG_PROFILER_SHINY
	UG_LOG("---------- Shiny Profiler Backtrace: ----------\n");
	Shiny::ProfileNode *p = Shiny::ProfileManager::instance._curNode;
	size_t i=0;

	while(p != &Shiny::ProfileManager::instance.rootNode)
	{
		const char *name = p->zone->name;
		if(name[0] == '@') name = "LUA Script";
		UG_LOG(std::setw(3) << i++ << std::setw(50) << name << "\t" << p->zone->file << " :" << p->zone->line << "\n");
		p = p->parent;
	}
	UG_LOG("\n");
#endif
#endif
}



string get_gcc_backtrace()
{
#ifdef UG_POSIX
	stringstream ss;
	void *array[100];
	size_t size;
	char **strings;

	size = backtrace (array, 100);
	strings = backtrace_symbols (array, size);

	for (size_t i = 0; i < size; i++)
		ss << i << ":\n" << ug::demangle_block(strings[i]);
	free (strings);
	ss << "\n";
	return ss.str();
#else
	return "";
#endif	// UG_POSIX
}

void gcc_backtrace()
{
#ifdef UG_POSIX
	UG_LOG("--------------- GCC Backtrace: ----------------\n");
	UG_LOG(get_gcc_backtrace());
#endif	// UG_POSIX
}



/*
 * This function is meant to put out some more meaningful debug data.
 * It prints the call backtrace and Profiler backtrace if available.
 * It also demangles the C++ function names so you get more a clue
 * what's going on.
 */
void ug_backtrace()
{
	lua_backtrace();
	shiny_backtrace();
	gcc_backtrace();

}

/// put a breakpoint here to break on UG_ASSERT or UG_THROW
/// b ug_assert_or_error
void ug_assert_or_error()
{
	// intentionally left empty
}


/// called whenever UG_ASSERT is called.
void ug_assert_failed()
{
	ug_backtrace();
	ug_assert_or_error();
#ifdef UG_PARALLEL
	pcl::Abort();
#endif
}

/// called whenever UG_THROW or UG_THROW_REGISTRY_ERROR is called.
void ug_throw_error()
{
	//ug_backtrace();
	ug_assert_or_error();
}
