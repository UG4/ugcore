/*
 * log_impl.h
 *
 *  Created on: 10.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__COMMON__LOG_IMPL__
#define __H__COMMON__LOG_IMPL__


namespace ug{

inline
std::ostream&
LogAssistant::
debug_logger()
{
	#ifdef UG_LOG_TO_FILE
	return file_logger();
	#else
	return std::cout;
	#endif
}

inline
std::ostream&
LogAssistant::
logger()
{
	#ifdef UG_LOG_TO_FILE
	return file_logger();
	#else
	return std::cout;
	#endif
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
	static LogAssistant log;
	return log;
}

} // end namespace ug

#endif /* __H__COMMON__LOG_IMPL__ */
