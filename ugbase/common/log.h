/*
 * log.h
 *
 *  Created on: 10.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__COMMON__LOG__
#define __H__COMMON__LOG__

#include <iostream>
#include <fstream>

namespace ug{

/* LogAssistant
 *
 * This class provides infrastructure for logging. It separates the log messages in different levels:
 *
 * logger() returns the output stream for normal User information output.
 * debug_logger() returns the output stream for debug informations.
 *
 * Debug messages are grouped by tags and debug levels. It is intended, that only messages are printed,
 * when the current level is equal or greater than the debug level chosen in the code.
 */
class LogAssistant
{
	public:
		// different tags to distingish several parts of the program, that
		// may be debugged separately
		enum Tags
		{
			MAIN = 0,
			LIB_GRID,
			LIB_GRID_REFINER,
			LIB_DISC,
			LIB_ALG,
			NUM_TAGS
		};

	public:
		// Constructor, resets all debug levels to -1
		LogAssistant();

		// returns the debug output stream
		inline std::ostream& debug_logger();

		// returns the normal output stream
		inline std::ostream& logger();

		// sets the debug level of all tags to 'lev'
		bool set_debug_levels(int lev);

		// sets the debug level of Tag 'tag' to level 'lev'
		bool set_debug_level(Tags tags, int lev);

		// returns the debug level of Tag 'tag'
		inline int get_debug_level(Tags tag);

	protected:
		// returns a file output stream for file logging (is used by logger(), debug_logger())
		std::ostream& file_logger();

	protected:
		// debug levels of tags
		int m_TagLevel[NUM_TAGS];
};

// returns singleton instance of LogAssistant
inline LogAssistant& GetLogAssistant();

} // end namespace ug


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// DEBUG LOG
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// Usage:
/* The following macros can be used to control debug messages.
 * To use them the define 'UG_ENABLE_DEBUG_LOGS' must be set. Otherwise nothing will be done (no runtime overhead).
 *
 *
 * UG_SET_DEBUG_LEVEL(tag, level)  	- sets the debug levels of Tag 'tag' to level 'level'
 * UG_RESET_DEBUG_LEVELS()			- sets the debug level of all Tags to -1
 * UG_SET_DEBUG_LEVELS(level)		- sets the debug levels of all Tags to level 'level'
 * UG_DLOG(tag, level, msg)			- prints the message "msg" to the debug log stream, if current level of Tag 'tag' is >= level
 *
 *	Example:
 *
	UG_RESET_DEBUG_LEVELS();
	UG_DLOG(MAIN, 0, "DLOG on level 0.\n");		// no message printed
	UG_SET_DEBUG_LEVEL(MAIN, 0);
	UG_DLOG(MAIN, 0, "DLOG on lebel 0.\n");		// message printed
	UG_DLOG(MAIN, 1, "DLOG on level 1.\n");		// no message printed
 */


#ifdef UG_ENABLE_DEBUG_LOGS
	#define UG_SET_DEBUG_LEVEL(tag, level)		{ug::GetLogAssistant().set_debug_level(ug::LogAssistant::tag, level);}
	#define UG_RESET_DEBUG_LEVELS()				{ug::GetLogAssistant().set_debug_levels(-1);}
	#define UG_SET_DEBUG_LEVELS(level)			{ug::GetLogAssistant().set_debug_levels(level);}
	#define UG_DLOG(tag, level, msg)			{if(ug::GetLogAssistant().get_debug_level(ug::LogAssistant::tag) >= level) {ug::GetLogAssistant().debug_logger() << msg; ug::GetLogAssistant().debug_logger().flush();}}
#else
	#define UG_SET_DEBUG_LEVEL(tag, level)		{}
	#define UG_RESET_DEBUG_LEVELS()				{}
	#define UG_SET_DEBUG_LEVELS(level)			{}
	#define UG_DLOG(tag, level, msg)			{}
#endif

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// WARNING LOG
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// Usage:
/* The following macro can be used to print warning messages.
 * To use it the define 'UG_ENABLE_WARNINGS' must be set. Otherwise nothing will be done (no runtime overhead).
 *
 * UG_WARNING(msg)  				- prints a warning to the normal output stream
 */


#ifdef UG_ENABLE_WARNINGS
	#define UG_WARNING(msg) {ug::GetLogAssistant().logger() << "UG_WARNING in " << __FILE__ << " at line " << __LINE__ << ": " << msg; ug::GetLogAssistant().logger().flush();}
#else
	#define UG_WARNING(msg) {}
#endif


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// LOG
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// Usage:
/*
 * UG_LOG(msg)  					- prints a message to the normal output stream
 */

#define UG_LOG(msg) {ug::GetLogAssistant().logger() << msg; ug::GetLogAssistant().logger().flush();}

// include implementation
#include "log_impl.h"

#endif /* __H__COMMON__LOG__ */
