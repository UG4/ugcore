/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

/*
 * This file is compiled if
 * cmake -DPROFILE_MEMORY=ON ..
 */


#include "memtracker.h"
#include "common/log.h"
#include "assert.h"
#include "common/util/string_util.h"
#include <map>


using namespace std;

namespace ug{

/////////////////////////////////////////////////////////////////////////////

map<const Shiny::ProfileNode *, size_t> selfmem;
map<const Shiny::ProfileNode *, size_t> totalmem;

class MemTrackerStruct
{
public:
	size_t size;
	Shiny::ProfileNode *pn;
	MemTrackerStruct()
	{

	}
	MemTrackerStruct(int s)
	{
		size = s;
		pn = Shiny::ProfileManager::instance._curNode;
		selfmem[pn]+=size;
	}
};
typedef map<void *, MemTrackerStruct> MemTrackerMap;


bool bMemTracker=false;
size_t allocated = 0;
MemTrackerMap memTracker;

/////////////////////////////////////////////////////////////////////////////

void *get_mem(size_t size)
{
	void *p = malloc( size );
	if(bMemTracker)
	{
		bMemTracker = false;
		MemTrackerStruct s(size);
		memTracker[p] = s;
		allocated+=size;
		bMemTracker = true;
	}
	if(p == NULL) throw std::bad_alloc();
	return p;
}


void release_mem(void *p)
{
	free(p);
	if(bMemTracker)
	{
		bMemTracker=false;
		MemTrackerMap::iterator it = memTracker.find(p);
		memTracker.erase(it);
		bMemTracker=true;
	}
}

} // namespace ug


void operator delete( void *p ) throw()
{
	ug::release_mem(p);
}
void operator delete[]( void *p ) throw()
{
	ug::release_mem(p);
}
void* operator new[]( size_t size ) throw(std::bad_alloc)
{
	return ug::get_mem(size);
}
void* operator new( size_t size ) throw(std::bad_alloc)
{
	return ug::get_mem(size);
}

namespace ug{

void DisplayVacantMemory()
{
	bool b = EnableMemTracker(false);
	bMemTracker = false;
	for(MemTrackerMap::iterator it = memTracker.begin(); it != memTracker.end(); ++it)
	{
		MemTrackerStruct &s = (*it).second;
		UG_LOG("vacant memory: size = " << GetBytesSizeString(s.size) << ", file " << s.pn->zone->file << " : " << s.pn->zone->line << "\n");
	}
	EnableMemTracker(b);
}


void CalcTotalMem(const Shiny::ProfileNode *p)
{
	size_t total = selfmem[p];
	for(const Shiny::ProfileNode *c=p->firstChild; c != NULL; c=c->nextSibling)
	{
		if(totalmem.find(c) == totalmem.end())
			CalcTotalMem(c);
		total += totalmem[c];
		if(c==p->lastChild)
			break;
	}
	totalmem[p] = total;
}

void UpdateTotalMem()
{
	const Shiny::ProfileNode *node = &Shiny::ProfileManager::instance.rootNode;
	bool b = EnableMemTracker(false);
	totalmem.clear();
	CalcTotalMem(node);
	EnableMemTracker(b);
}

size_t GetSelfMem(const Shiny::ProfileNode *p)
{
	return selfmem[p];
}
size_t GetTotalMem(const Shiny::ProfileNode *p)
{
	return totalmem[p];
}

bool EnableMemTracker(bool b)
{
	bool b2 = bMemTracker;
	bMemTracker =b;
	return b2;
}
bool IsMemTrackerEnabled()
{
	return bMemTracker;
}

bool HasMemTracking()
{
	return true;
}

} // namespace ug



