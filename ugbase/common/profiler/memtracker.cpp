/*
 * memtracker.cpp
 *
 *  Created on: 30.04.2013
 *      Author: mrupp
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



