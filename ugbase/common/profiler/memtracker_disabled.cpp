/*
 * This file is compiled if
 * cmake -DPROFILE_MEMORY=ON ..
 */

#include "memtracker.h"
#include "common/log.h"
#include "assert.h"
#include <map>


using namespace std;

namespace ug{

#ifdef UG_PROFILER_SHINY
size_t GetSelfMem(const Shiny::ProfileNode *p)
{
	return 0;
}
size_t  GetTotalMem(const Shiny::ProfileNode *p)
{
	return 0;
}
#endif
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



