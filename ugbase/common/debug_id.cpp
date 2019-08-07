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

#include "debug_id.h"
#include "common/log.h"
#include "common/error.h"
#include "common/assert.h"
#include <iostream>
#include <string.h>
#include "util/string_util.h"

namespace ug{

/**
 * register the debug id.
 * NOTE: this function is called in the initialization of variables
 * of type DebugID, which are mostly global variables.
 * Be absolutely sure we are safe here, i.e.
 * we are not using other things which might not be initialized (like Log).
 */
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
	int slen = strlen(debugID);
	if(slen<=0) return false;
	if(debugID[slen-1] == '*')
	{
		for(size_t i=0; i<m_dbgLevelIdentifiers.size(); i++)
		{
			const char *name = m_dbgLevelIdentifiers[i].c_str();
			if(WildcardMatch(name, debugID))
			{
				set_debug_level(crc32(name), level);
				//UG_LOGN(name);
			}
		}
	}
	else if(set_debug_level(crc32(debugID), level) == false)
	{
		UG_LOG("DebugID " << debugID << " not registered.\n");
		return false;
	}
	return true;
}

/**
 * register the debug id.
 * NOTE: this function is called in the initialization of global variables
 * of type DebugID. Be absolutely sure we are safe here, i.e.
 * we are not using other things which might not be initialized (like Log).
 */
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
		// --> it should be, if we define it in this header!
		std::cout << "FATAL ERROR: DebugID "<<debugID<<" already registered." << std::endl;
		// note that this could be caused by double-registering libraries.
		UG_THROW("FATAL ERROR: DebugID "<<debugID<<" already registered.");
		return false;
	}
}

}
