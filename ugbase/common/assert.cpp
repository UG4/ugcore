/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include "common/log.h"
#include "common/profiler/profiler.h"
#include "assert.h"
#include <stdlib.h>

#ifdef UG_POSIX
#include "util/demangle.h"
#ifndef ANDROID
#include <execinfo.h>
#endif // ANDROID
#endif // UG_POSIX

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
#if defined UG_POSIX && !defined ANDROID
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
