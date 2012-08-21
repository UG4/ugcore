/*
 * assert.cpp
 *
 *  Created on: 05.07.2012
 *      Author: Martin Rupp
 */

#include "common/log.h"
#include "common/profiler/profiler.h"
#include <stdlib.h>
#ifdef UG_POSIX
#include <execinfo.h>
#include <cxxabi.h>
#endif
#include <string.h>

/*
 * This function is meant to put out some more meaningful assert data.
 * It demangles the C++ function names so you get more a clue
 * what's going on.
 */
void ug_assert_failed()
{
#ifdef UG_POSIX

#ifdef UG_PROFILER
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
	{
		UG_LOG(i << ":\n");
		int status;
		// library(mangledFunction+adress) [address]
		char *f1=strchr(strings[i], '(');
		char *f2=strchr(strings[i], '+');
		if(f1==NULL || f2 == NULL || f1 > f2)
		{ 	UG_LOG(strings[i] << "\n"); }
		else
		{
			*f1=0x00; *f2=0x00;
			char *realname = abi::__cxa_demangle(f1+1, 0, 0, &status);
			if(status==0)
			{   UG_LOG(strings[i] << "(" << realname << "+" << f2+1 << "\n");}
			else
			{	UG_LOG(strings[i] << "(" << f1+1 << "+" << f2+1 << "\n"); }
			free(realname);
		}
	}
	free (strings);
	UG_LOG("\n");

#endif	// UG_POSIX

}
