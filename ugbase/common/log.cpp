/*
 * log.cpp
 *
 *  Created on: 10.03.2010
 *      Author: andreasvogel
 */

#include "common/log.h"

namespace ug{

LogAssistant::
LogAssistant()
{
	set_debug_levels(-1);
}

bool
LogAssistant::
set_debug_levels(int lev)
{
	for(int i = 0; i < NUM_TAGS; ++i)
	{
		m_TagLevel[i] = lev;
	}
	return true;
}

bool
LogAssistant::
set_debug_level(Tags tags, int lev)
{
	m_TagLevel[tags] = lev;
	return true;
}

std::ostream&
LogAssistant::
file_logger()
{
	static std::ofstream out("uglog.log");
	return out;
}

}
