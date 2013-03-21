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

/*
 * This function is meant to put out some more meaningful debug data.
 * It prints the call backtrace and Profiler backtrace if available.
 * It also demangles the C++ function names so you get more a clue
 * what's going on.
 */
void ug_backtrace()
{
	
#ifdef UG_POSIX

#ifdef UG_PROFILER_SHINY
	UG_LOG("Profiler Backtrace:\n")
	Shiny::ProfileNode *p = Shiny::ProfileManager::instance._curNode;
	while(p != &Shiny::ProfileManager::instance.rootNode)
	{
		UG_LOG(p->zone->name << "\n");
		p = p->parent;
	}
	UG_LOG("\n");
#endif
	UG_LOG("GCC Backtrace:\n")
	void *array[100];
	size_t size;
	char **strings;
	size_t i;

	size = backtrace (array, 100);
	strings = backtrace_symbols (array, size);

	for (i = 0; i < size; i++)
		UG_LOG(i << ":\n" << demangle(strings[i]));
	free (strings);
	UG_LOG("\n");

#endif	// UG_POSIX

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
