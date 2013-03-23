/* 
 * File:   DebugID.cpp
 * Author: mrupp
 * 
 * Created on 22. Oktober 2012, 13:55
 */

#include "debug_id.h"
#include "common/log.h"
#include "common/error.h"
#include "common/assert.h"

namespace ug{

DebugID::DebugID(const char *str)
{
	m_hash = crc32(str);
	GetDebugIDManager().register_debug_id(str);
}

DebugIDManager& DebugIDManager::instance()
{
	static DebugIDManager m;
	return m;
}


bool DebugIDManager::
set_debug_levels(int lev)
{
	for(std::map<uint32, int>::iterator it=m_dbgLevels.begin(); it != m_dbgLevels.end(); ++it)
		(*it).second = lev;
	return true;
}

bool DebugIDManager::
set_debug_level(const char *debugID, int level)
{
	if(set_debug_level(crc32(debugID), level) == false)
	{
		UG_LOG("DebugID " << debugID << " not registered.\n");
		return false;
	}
	else return true;
}

bool DebugIDManager::
register_debug_id(const char *debugID)
{
	if(debug_id_registered(debugID) == false)
	{
		m_dbgLevelIdentifiers.push_back(std::string(debugID));
		m_dbgLevels[crc32(debugID)] = -1;
		return true;
	}
	else
	{
		// not quite clear if cout is defined yet.
		// note that this could be caused by double-registering libraries.
		UG_THROW("FATAL ERROR: DebugID "<<debugID<<" already registered.");
		return false;
	}
}

}
