/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Sebastian Reiter
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

#ifndef __H__UG__COMMON__LOG__
#define __H__UG__COMMON__LOG__

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "util/ostream_buffer_splitter.h"
#include "types.h"
#include "ug_config.h"

#include <vector>
#include <map>
#include "common/util/crc32.h"
#include "debug_id.h"


//	in order to support VRL logs, we're including bindings_vrl.h
//  this is necessary to get access to the JVM environment
#ifdef UG_FOR_VRL
		#include <sstream>
		#include "bindings/vrl/messaging.h"
#endif

namespace ug{

/// \addtogroup ugbase_common
/// \{

	const uint64 UNIT_KILO     = 1024;                   // 2^{10}
	const uint64 UNIT_MEGA     = UNIT_KILO * 1024;       // 2^{20} =                 1'048'576
	const uint64 UNIT_GIGA     = UNIT_MEGA * 1024;       // 2^{30} =             1'073'741'824
	const uint64 UNIT_TERA     = UNIT_GIGA * 1024ll;     // 2^{40} =         1'099'511'627'776
	const uint64 UNIT_PETA     = UNIT_TERA * 1024ll;     // 2^{50} =     1'125'899'906'842'624
	const uint64 UNIT_EXA      = UNIT_PETA * 1024ll;     // 2^{60} = 1'152'921'504'606'846'976
	//const uint64_t UNIT_ZETTA    = UNIT_EXA  * 1024ll;   // 2^{70} -- too big for 64 bit 'long long int'!
	
	const uint64 UNIT_KILO_SI  = 1000;                   // 10^{ 3}
	const uint64 UNIT_MEGA_SI  = UNIT_KILO_SI * 1000;    // 10^{ 6}
	const uint64 UNIT_GIGA_SI  = UNIT_MEGA_SI * 1000;    // 10^{ 9}
	const uint64 UNIT_TERA_SI  = UNIT_GIGA_SI * 1000ll;  // 10^{12}
	const uint64 UNIT_PETA_SI  = UNIT_TERA_SI * 1000ll;  // 10^{15}
	const uint64 UNIT_EXA_SI   = UNIT_PETA_SI * 1000ll;  // 10^{18}
	//const uint64_t UNIT_ZETTA_SI = UNIT_EXA_SI  * 1000ll;// 10^{21} -- too big for 64 bit 'long long int'!



// LogAssistant
/**
 * This class provides infrastructure for logging. It separates the log messages
 * in different levels:
 *
 * logger() returns the output stream for normal User information output.
 * debug_logger() returns the output stream for debug informations.
 *
 * Debug messages are grouped by tags and debug levels. It is intended, that
 * only messages are printed, when the current level is equal or greater than
 * the debug level chosen in the code.
 *
 * Please note that this class operates on std::clog. Thus, if output-options
 * are changed (e.g. file-logging enabled) the stream buffer on which clog
 * operates will change too.
 */
class UG_API LogAssistant
{
	public:
	///	returns a reference to the single instance of LogAssistant
		static LogAssistant& instance();

	///	flushes all buffers before deletion
		~LogAssistant();

	/// enables or disables file output.
	/** a filename can be specified. Default is 'uglog.log'.
	 *	Please note that only the filename given at the first call is considered.
	 *	Filelogging is disabled by default.
	 */
		bool enable_file_output(bool bEnable, const char* filename = "uglog.log");

	/// renames output file if opened.
		bool rename_log_file(const char * newname);

	///	enables or disables terminal output.
	/**	terminal output is enbled by default.*/
		bool enable_terminal_output(bool bEnable);

	/// returns the debug output stream
		inline std::ostream& debug_logger();

	/// returns the normal output stream
		inline std::ostream& logger();

	///	returns the error output stream
	/**	Note that error-logs are not immediately visible. Instead they are
	 * gathered and can be printed in one block using LogAssistant::flush_error_log().*/
		inline std::ostream& error_logger();

	///	outputs all gathered error-logs to the standard logger and clears the error log.
		void flush_error_log();

	/// sets the debug level of all tags to 'lev'
		bool set_debug_levels(int lev)
		{
			return GetDebugIDManager().set_debug_levels(lev);
		}

	/// returns the debug level of debugID
		inline int get_debug_level(DebugID &debugID) const
		{
			return GetDebugIDManager().get_debug_level(debugID);
		}
	/// returns the debug level of debugID
		inline int get_debug_level(const char *debugID) const
		{
			return GetDebugIDManager().get_debug_level(debugID);
		}
		/// returns the debug level of debugID
		int get_debug_level_noninline(const char *debugID) const;
		
		
	/// sets the debug level of debugID
		inline bool set_debug_level(DebugID &debugID, int level)
		{
			return GetDebugIDManager().set_debug_level(debugID, level);
		}
		
	/// returns the debug level of debugID
		inline bool set_debug_level(const char* debugID, int level)
		{
			return GetDebugIDManager().set_debug_level(debugID, level);
		}
		
		bool set_debug_level_noninline(const char* debugID, int level);
		
	/// returns the debug level of debugID
		inline std::string get_registered_debug_IDs()
		{
			return GetDebugIDManager().get_registered_debug_IDs();
		}

	///	sets the output-process in a parallel environment. Default is 0.
	/**	pass a procRank of -1 to enable output on all processes.*/
		void set_output_process(int procRank);

	///	returns the rank of the current output-process or -1 if all procs perform output.
		inline int get_output_process();

	///	returns true if the current process is an output process.
	/**	This is always true, if the application is executed in a serial environment.*/
		bool is_output_process();

	///	returns the process rank of the underlying process (same as pcl::ProcRank)
		int get_process_rank();

		void flush();

	protected:
	///	updates and sets stream buffers based on current options.
	/**	Note that this method changes the buffer on which clog works.*/
		void update_ostream();

	///	opens the local logFile, if the process is the output-process
		bool open_logfile();

	private:
	/// Constructor, resets all debug levels to -1
	/**	Constructor is private since only one instance of LogAssistant may exist.*/
		LogAssistant();

	///	Prevent copying through private copy constructor.
		LogAssistant(const LogAssistant&);

	///	Performs some initialization
		void init();
	private:
	//	streams
		std::streambuf*	m_emptyBuf;
		std::streambuf*	m_logBuf;
		std::streambuf*	m_fileBuf;
		std::streambuf*	m_splitBuf;

		EmptyStreamBuffer		m_emptyBufInst;
		OStreamBufferSplitter	m_splitBufInst;

		std::string 			m_logFileName;
		std::ofstream			m_fileStream;
		std::stringstream		m_errStream;

		bool m_terminalOutputEnabled;
		bool m_fileOutputEnabled;

		int m_outputProc;

};

// returns singleton instance of LogAssistant
inline LogAssistant& GetLogAssistant();

/// returns number 'size' in a more human readable format (using IEC binary prefixes)
inline std::string ConvertNumber (uint64_t size, unsigned int width,
                                  unsigned int numDisplayedDigits);

/// returns number 'size' in a more human readable format (using SI prefixes)
inline std::string ConvertNumberSI (uint64_t size, unsigned int width,
                                    unsigned int numDisplayedDigits);

// end group ugbase_common
/// \}

} // end namespace ug


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// DEBUG LOG
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/// \addtogroup ugbase_common
/// \{

// Usage:
/* The following macros can be used to control debug messages.
 * To use them the define 'UG_ENABLE_DEBUG_LOGS' must be set. Otherwise nothing will be done (no runtime overhead).
 *
 * UG_SET_DEBUG_LEVEL(__debugID__, level)  	- sets the debug levels of DebugID __debugID__ to level 'level'
 * UG_RESET_DEBUG_LEVELS()					- sets the debug level of all Tags to -1
 * UG_SET_DEBUG_LEVELS(level)				- sets the debug levels of all DebugIDs to level 'level'
 * UG_DLOG(__debugID__, level, msg)			- prints the message "msg" to the debug log stream, if current level of DebugID __debugID__ is >= level
 *
 * __debugID__ has to be of type DebugID.
 * For a list of standard DebugIDs and how to define your own, see debug_id.h .
 *
 *	Example:
 *
	UG_RESET_DEBUG_LEVELS();
	UG_DLOG(MAIN, 0, "DLOG on level 0.\n");		// no message printed
	UG_SET_DEBUG_LEVEL(MAIN, 0);
	UG_DLOG(MAIN, 0, "DLOG on level 0.\n");		// message printed
	UG_DLOG(MAIN, 1, "DLOG on level 1.\n");		// no message printed
 */
#ifdef UG_ENABLE_DEBUG_LOGS
	#define UG_SET_DEBUG_LEVEL(__debugID__, level)		{ug::GetDebugIDManager().set_debug_level(__debugID__, level);}
	#define UG_RESET_DEBUG_LEVELS()						{ug::GetDebugIDManager().set_debug_levels(-1);}
	#define UG_SET_DEBUG_LEVELS(level)					{ug::GetDebugIDManager().set_debug_levels(level);}
	#define UG_DEBUG_BEGIN(__debugID__, level)			{ if(ug::GetDebugIDManager().get_debug_level(__debugID__) >= level) {
	#define UG_DEBUG_END(__debugID__, level)			}; }
	#define IF_DEBUG(__debugID__, level) 				if(ug::GetDebugIDManager().get_debug_level(__debugID__) >= level)

	#define UG_DLOG(__debugID__, level, msg) \
			{if(ug::GetDebugIDManager().get_debug_level(__debugID__) >= level)\
			{ug::GetLogAssistant().debug_logger() << msg; ug::GetLogAssistant().debug_logger().flush();}}

	#define UG_DLOGN(__debugID__, level, msg)			UG_DLOG(__debugID__, level, msg << "\n");

	#define UG_DLOG_ALL_PROCS(__debugID__, level, msg)	\
			{if(ug::GetDebugIDManager().get_debug_level(__debugID__) >= level)\
			{ug::LogAssistant& la = ug::GetLogAssistant(); int op = la.get_output_process();\
			la.set_output_process(-1); la.debug_logger() << "[Proc " << la.get_process_rank() << "]: "\
			<< msg << std::flush; la.set_output_process(op);}}
#else
	#define UG_SET_DEBUG_LEVEL(__debugID__, level)		{}
	#define UG_RESET_DEBUG_LEVELS()						{}
	#define UG_SET_DEBUG_LEVELS(level)					{}
	#define UG_DLOG(__debugID__, level, msg)			{}
	#define UG_DLOGN(__debugID__, level, msg)			{}
	#define UG_DLOG_ALL_PROCS(__debugID__, level, msg)	{}
	#define UG_DEBUG_BEGIN(__debugID__, level)			{ if(1==0) {
	#define UG_DEBUG_END(__debugID__, level)			}; }
	#define IF_DEBUG(__debugID__, level)				if(1==0)
#endif

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// WARNING LOG
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Warning and debug logs are disabled when compiling for release
// defines (currently here, should be compiling options)
#ifndef NDEBUG
	#define UG_ENABLE_WARNINGS
#else
	#undef UG_ENABLE_WARNINGS
#endif

// Usage:
/* The following macro can be used to print warning messages.
 * To use it the define 'UG_ENABLE_WARNINGS' must be set. Otherwise nothing will
 * be done (no runtime overhead).
 *
 * UG_WARNING(msg)  			- prints a warning to the normal output stream
 */
#ifdef UG_ENABLE_WARNINGS
	#define UG_WARNING(msg) {ug::GetLogAssistant().logger() << "UG_WARNING in "\
								<< __FILE__ << " at line " << __LINE__ << ": " \
								<< msg << std::flush;}
#else
	#define UG_WARNING(msg) {}
#endif


// Usage:
/* The following macro can be used to print warning messages if condition is met.
 * To use it the define 'UG_ENABLE_WARNINGS' must be set. Otherwise nothing will
 * be done (no runtime overhead).
 *
 * UG_WARNING(cond, msg)  		- prints a warning to the normal output stream
 */
#ifdef UG_ENABLE_WARNINGS
	#define UG_COND_WARNING(cond, msg) { if (cond) {ug::GetLogAssistant().logger() << "UG_WARNING in "\
							<< __FILE__ << " at line " << __LINE__ << ": " \
							<< msg << std::flush;} }
#else
	#define UG_COND_WARNING(cond, msg) {}
#endif





////////////////////////////////////////////////////////////////////////////////
// LOG
/**
 * UG_LOG(msg)  		- prints a message to the normal output stream
 */
#ifdef UG_FOR_VRL
	#define VRL_LOG(msg) {std::stringstream __ss; __ss << "" << msg;\
						  ug::vrl::MessageBuffer::addMessage(__ss.str());}
#else
	#define VRL_LOG(msg)
#endif

#define UG_LOG(msg) {ug::GetLogAssistant().logger() << msg << std::flush; VRL_LOG(msg);}
#define UG_COND_LOG(cond, msg) { if (cond) { UG_LOG(msg); } }
#define UG_LOGN(msg) UG_LOG(msg << "\n")
#define UG_COND_LOGN(cond, msg) { if (cond) { UG_LOGN(msg); } }
#define UG_LOG_ALL_PROCS(msg) {ug::LogAssistant& la = ug::GetLogAssistant();\
							   int op = la.get_output_process(); la.set_output_process(-1);\
							   la.logger() << "[Proc " << std::setw(3) << la.get_process_rank() << "]: "\
							   << msg << std::flush; VRL_LOG(msg); la.set_output_process(op);}

////////////////////////////////////////////////////////////////////////////////
// LOG
/**
 * UG_ERR_LOG(msg)  prints msg to LogAssistant::error_logger().
 * 					Those messages are not immediately flushed. Instead they can
 * 					be flushed to the standard logger using LogAssistant::flush_error_log
 */
#define UG_ERR_LOG(msg) {ug::GetLogAssistant().error_logger() << msg; VRL_LOG(msg);}

// include implementation
#include "log_impl.h"

// end group ugbase_common
/// \}

#endif