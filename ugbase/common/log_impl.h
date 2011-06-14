/*
 * log_impl.h
 *
 *  Created on: 10.03.2010
 *      Author: andreasvogel, sebastianreiter
 */

#ifndef __H__COMMON__LOG_IMPL__
#define __H__COMMON__LOG_IMPL__

#include <iostream>

namespace ug{

inline
std::ostream&
LogAssistant::
debug_logger()
{
	return std::cout;
}

inline
std::ostream&
LogAssistant::
logger()
{
	return std::cout;
}

inline
int
LogAssistant::
get_debug_level(Tags tag)
{
	return m_TagLevel[tag];
}


inline
LogAssistant&
GetLogAssistant()
{
	return LogAssistant::instance();
}

} // end namespace ug

#endif /* __H__COMMON__LOG_IMPL__ */
