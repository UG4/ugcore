/*
 * log.h
 *
 *  Created on: 10.03.2010
 *      Author: andreasvogel, sebastianreiter
 */

#ifndef __H__COMMON__LOG__
#define __H__COMMON__LOG__

#include <iostream>
#include <fstream>
#include "util/ostream_util.h"

//	in order to support parallel logs, we're including pcl.h
//	you can set the output process using pcl::SetOutputProcRank(int rank).
#ifdef UG_PARALLEL
	#include "pcl/pcl_base.h"
#endif

//	in order to support VRL logs, we're including bindings_vrl.h
//  this is necessary to get access to the JVM environment
#ifdef FOR_VRL
		#include <sstream>
		#include "bindings_vrl/messaging.h"
#endif

namespace ug{
// LogAssistant
/**
 * This class provides infrastructure for logging. It separates the log messages in different levels:
 *
 * logger() returns the output stream for normal User information output.
 * debug_logger() returns the output stream for debug informations.
 *
 * Debug messages are grouped by tags and debug levels. It is intended, that only messages are printed,
 * when the current level is equal or greater than the debug level chosen in the code.
 *
 * Please note that this class operates on std::clog. Thus, if output-options are changed
 * (e.g. file-logging enabled) the stream buffer on which clog operates will change too.
 */
class LogAssistant
{
	public:
		// different tags to distingish several parts of the program, that
		// may be debugged separately
		enum Tags
		{
			MAIN = 0,
			APP,
			LIB_GRID,
			LIB_GRID_REFINER,
			LIB_DISC,
			LIB_DISC_ASSEMBLE,
			LIB_DISC_D3F,
			LIB_DISC_MULTIGRID,
			LIB_DISC_NEWTON,
			LIB_DISC_LINKER,
			LIB_DISC_TRANSFER,
			LIB_DISC_DISCRETE_FUNCTION,
			LIB_DISC_OUTPUT,
			LIB_DISC_OPERATOR_INVERSE,
			LIB_ALG_LINEAR_OPERATOR,
			LIB_ALG_LINEAR_SOLVER,
			LIB_ALG_VECTOR,
			LIB_ALG_MATRIX,
			NUM_TAGS
		};

	public:
	///	returns a reference to the single instance of LogAssistant
		static LogAssistant& instance();

	/// enables or disables file output.
	/** a filename can be specified. Default is 'uglog.log'.
	 *	Please note that only the filename given at the first call is considered.
	 *	Filelogging is disabled by default.
	 */
		bool enable_file_output(bool bEnable, const char* filename = "uglog.log");

	///	enables or disables terminal output.
	/**	terminal output is enbled by default.*/
		bool enable_terminal_output(bool bEnable);

	/// returns the debug output stream
		inline std::ostream& debug_logger();

	/// returns the normal output stream
		inline std::ostream& logger();

	/// sets the debug level of all tags to 'lev'
		bool set_debug_levels(int lev);

	/// sets the debug level of Tag 'tag' to level 'lev'
		bool set_debug_level(Tags tags, int lev);

	/// returns the debug level of Tag 'tag'
		inline int get_debug_level(Tags tag);

	protected:
	///	updates and sets stream buffers based on current options.
	/**	Note that this method changes the buffer on which clog works.*/
		void update_ostream();

	private:
	/// Constructor, resets all debug levels to -1
	/**	Constructor is private since only one instance of LogAssistant may exist.*/
		LogAssistant();

	///	Prevent copying through private copy constructor.
		LogAssistant(const LogAssistant& la)	{}

	private:
	//	streams
		std::streambuf*	m_emptyBuf;
		std::streambuf*	m_logBuf;
		std::streambuf*	m_fileBuf;
		std::streambuf*	m_splitBuf;

		OStreamBufferEmpty		m_emptyBufInst;
		OStreamBufferSplitter	m_splitBufInst;

		std::ofstream			m_fileStream;

		bool m_terminalOutputEnabled;
		bool m_fileOutputEnabled;

	/// debug levels of tags
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
	#define UG_DEBUG_BEGIN(tag, level)			{ if(ug::GetLogAssistant().get_debug_level(ug::LogAssistant::tag) >= level) {
	#define UG_DEBUG_END(tag, level)			}; }

	#ifdef UG_PARALLEL
		#define UG_DLOG(tag, level, msg)			{if(pcl::IsOutputProc()){\
														if(ug::GetLogAssistant().get_debug_level(ug::LogAssistant::tag) >= level)\
														{ug::GetLogAssistant().debug_logger() << msg; ug::GetLogAssistant().debug_logger().flush();}}}
		#define UG_DLOG_ALL_PROCS(tag, level, msg)	{if(ug::GetLogAssistant().get_debug_level(ug::LogAssistant::tag) >= level)\
													{ug::GetLogAssistant().debug_logger() << "[Proc " << pcl::GetProcRank() << "]: "; \
														ug::GetLogAssistant().debug_logger() << msg; ug::GetLogAssistant().debug_logger().flush();}}
	#else
		#define UG_DLOG(tag, level, msg)			{if(ug::GetLogAssistant().get_debug_level(ug::LogAssistant::tag) >= level)\
													 {ug::GetLogAssistant().debug_logger() << msg; ug::GetLogAssistant().debug_logger().flush();}}
		#define UG_DLOG_ALL_PROCS(tag, level, msg)	{if(ug::GetLogAssistant().get_debug_level(ug::LogAssistant::tag) >= level)\
													{ug::GetLogAssistant().debug_logger() << "[Proc 0]: "; \
													ug::GetLogAssistant().debug_logger() << msg; ug::GetLogAssistant().debug_logger().flush();}}
	#endif
#else
	#define UG_SET_DEBUG_LEVEL(tag, level)		{}
	#define UG_RESET_DEBUG_LEVELS()				{}
	#define UG_SET_DEBUG_LEVELS(level)			{}
	#define UG_DLOG(tag, level, msg)			{}
	#define UG_DLOG_ALL_PROCS(tag, level, msg)	{}
	#define UG_DEBUG_BEGIN(tag, level)			{ if(1==0) {
	#define UG_DEBUG_END(tag, level)			}; }
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
/**
 * UG_LOG(msg)  					- prints a message to the normal output stream
 * If ug is compiled in a parallel environment (UG_PARALLEL is defined),
 * UG_LOG will use PCLLOG to output its data.
 */
#ifdef FOR_VRL
	#define VRL_LOG(msg) {std::stringstream ss;ss << "<!--UG4-->" << msg;\
                          ug::vrl::MessageBuffer::addMessage(ss.str());}
#else
	#define VRL_LOG(msg)
#endif

#ifdef UG_PARALLEL
	#define UG_LOG(msg) {if(pcl::IsOutputProc())\
						{ug::GetLogAssistant().logger() << msg; VRL_LOG(msg);\
						 ug::GetLogAssistant().logger().flush();}}
	#define UG_LOG_ALL_PROCS(msg) {ug::GetLogAssistant().logger() << "[Proc " << pcl::GetProcRank() << "]: "; \
						           ug::GetLogAssistant().logger() << msg;\
						           VRL_LOG(msg);\
						           ug::GetLogAssistant().logger().flush();}
#else
	#define UG_LOG(msg) {ug::GetLogAssistant().logger() << msg; VRL_LOG(msg);\
						 ug::GetLogAssistant().logger().flush();}
	#define UG_LOG_ALL_PROCS(msg) {ug::GetLogAssistant().logger() << "[Proc 0]: "; \
						           ug::GetLogAssistant().logger() << msg;\
						           VRL_LOG(msg);\
						           ug::GetLogAssistant().logger().flush();}
#endif


//#define UG_LOG(msg) {std::stringstream ss;ss << msg; ug::vrl::soutPrintln(ss.str());}


// include implementation
#include "log_impl.h"

#endif /* __H__COMMON__LOG__ */
