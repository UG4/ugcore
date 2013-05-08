/*
 * memtracker_disabled.cpp
 *
 *  Created on: 08.05.2013
 *      Author: mrupp
 */

#include "memtracker.h"
#include "common/log.h"
#include "assert.h"
#include <map>


using namespace std;

namespace ug{

string GetBytesSize(size_t s, int length)
{
	stringstream ss;
	if(length!=0)
		ss << setw(length-3);
	if(s > 1024*1024*1024)
		ss << s/(1024*1024*1024.0) << " Gb";
	else if(s > 1024*1024)
		ss << s/(1024*1024.0) << " Mb";
	else if(s > 1024)
			ss << s/(1024.0) << " kb";
	else if(length == 0)
		ss << s << " b";
	else
		ss << s << " b ";
	return ss.str();
}

size_t GetSelfMem(const Shiny::ProfileNode *p)
{
	return 0;
}
size_t  GetTotalMem(const Shiny::ProfileNode *p)
{
	return 0;
}

bool EnableMemTracker(bool b)
{
	return false;
}
bool IsMemTrackerEnabled()
{
	return false;
}

bool HasMemTracking()
{
	return false;
}
void DisplayVacantMemory()
{
}
void UpdateTotalMem()
{

}

} // namespace ug



