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
#include <execinfo.h>
#include <cxxabi.h>
#endif
#include <string>
#include <sstream>
#ifdef UG_FOR_LUA
#include "bindings/lua/info_commands.h"
#include <map>
namespace ug{
namespace script{
bool IsLUADebug();
}
}
#endif

using namespace std;




#ifdef UG_POSIX
/**
 * demangles C++ function names like _ZZ12ug_backtracev = ug_backtrace().
 * also demangles them when they appear "in between", make sure they are
 * seperated by ' ', '\n' or '\t' and start with _.
 * @param str mangled string, containing stuff like _ZZ12ug_backtracev
 * @return the demangled string
 */
string demangle(const char *str)
{
	stringstream ss;
	char lastc=0x00;
    string s;
	int status;
	for(char c = *str; c != 0x00; c = *(++str))	
	{
		// mangled names start with _ . consider only if last sign was space or tab or newline
		if(c == '_' && (lastc==' ' || lastc == '\n' || lastc == '\t'))
		{
			s = "_";
			// add all characters to the string until space, tab or newline
			for(c = *(++str); c != 0x00; c = *(++str))			
			{
				if(c == ' ' || c == '\n' || c == '\t')
					break;
				s += c;
			}
			// some compilers add an additional _ in front. skip it.
			const char *p = s.c_str();
			if(s.length() > 2 && s[1] == '_') p = p+1;
			char *realname = abi::__cxa_demangle(p, 0, 0, &status);
			if(status==0)
				ss << realname; // demangle successfull
			else
				ss << s; // demangling failed, print normal string
			free(realname);
		}
		ss << c;
		lastc =c;
	}
	return ss.str();
}
#endif

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
		ug::bridge::LuaStackTrace();
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
		ss << i << ":\n" << demangle(strings[i]);
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
	// intentially left empty
}


/// called whenever UG_ASSERT is called.
void ug_assert_failed()
{
	ug_backtrace();
	ug_assert_or_error();
}

/// called whenever UG_THROW or UG_THROW_REGISTRY_ERROR is called.
void ug_throw_error()
{
	//ug_backtrace();
	ug_assert_or_error();
}
