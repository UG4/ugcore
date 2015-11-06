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

#include <iostream>
#include "common/log.h"
#include "common/profiler/profiler.h"
#include <cassert>

//	in order to support parallel logs, we're including pcl.h
#ifdef UG_PARALLEL
	#include "pcl/pcl_base.h"
#endif


using namespace std;

namespace ug{

// predefined standard DebugIDs. todo: move some of those in the appropriate modules,
// syntax for subgroups ?
DebugID MAIN("MAIN"),
		APP("APP"),
		LIB_GRID("LIB_GRID"),
		LIB_GRID_REFINER("LIB_GRID_REFINER"),
		LIB_DISC("LIB_DISC"),
		LIB_DISC_ASSEMBLE("LIB_DISC_ASSEMBLE"),
		LIB_DISC_D3F("LIB_DISC_D3F"),
		LIB_DISC_MULTIGRID("LIB_DISC_MULTIGRID"),
		LIB_DISC_NEWTON("LIB_DISC_NEWTON"),
		LIB_DISC_LINKER("LIB_DISC_LINKER"),
		LIB_DISC_TRANSFER("LIB_DISC_TRANSFER"),
		LIB_DISC_DISCRETE_FUNCTION("LIB_DISC_DISCRETE_FUNCTION"),
		LIB_DISC_OUTPUT("LIB_DISC_OUTPUT"),
		LIB_DISC_OPERATOR_INVERSE("LIB_DISC_OPERATOR_INVERSE"),
		LIB_ALG_LINEAR_OPERATOR("LIB_ALG_LINEAR_OPERATOR"),
		LIB_ALG_LINEAR_SOLVER("LIB_ALG_LINEAR_SOLVER"),
		LIB_ALG_VECTOR("LIB_ALG_VECTOR"),
		LIB_ALG_MATRIX("LIB_ALG_MATRIX"),
		LIB_ALG_AMG("LIB_ALG_AMG"),
		LIB_PCL("LIB_PCL");


///////////////////////////////////////////////////////////////////////////

LogAssistant::
LogAssistant() :
	m_terminalOutputEnabled(true),
	m_fileOutputEnabled(false),
	m_outputProc(0)
{
}

LogAssistant::LogAssistant(const LogAssistant&)
{
}


LogAssistant::~LogAssistant()
{
	flush();
	if(m_fileStream.is_open())
		m_fileStream.close();
}


void LogAssistant::init()
{
	PROFILE_FUNC();
//	originally this was placed in the constructor.
//	However, an internal compiler error in MinGW (windows)
//	appears then. Having moved the stuff here into this init
//	method seems to have helped.
	set_debug_levels(-1); // 0 or -1 ?
	m_emptyBuf = &m_emptyBufInst;
	m_splitBuf = &m_splitBufInst;
	m_logBuf = cout.rdbuf();
	m_fileBuf = m_fileStream.rdbuf();
	m_splitBufInst.set_buffers(m_logBuf, m_fileBuf);
	update_ostream();
}

LogAssistant& LogAssistant::
instance()
{
	static LogAssistant log;
	
//	This special initialization is performed to avoid an internal compiler error
//	with MinGW. I don't know what causes that error, however using a init method
//	instead of the constructor seems to help.
	static bool initialized = false;
	if(!initialized){
		initialized = true;
		log.init();
	}
	
	return log;
}

void LogAssistant::flush()
{
	logger().flush();
	m_splitBufInst.flush();
	if(m_fileStream.is_open())
		m_fileStream.flush();
	if(m_terminalOutputEnabled)
		fflush(stdout);
}


bool LogAssistant::
enable_file_output(bool bEnable, const char* filename)
{
	PROFILE_FUNC();
	m_logFileName = filename;
	if(bEnable){
		if(!open_logfile()){
			update_ostream();
			return false;
		}
	}
	else
		m_fileStream.close();

	m_fileOutputEnabled = bEnable;
	update_ostream();
	return true;
}

bool LogAssistant::
rename_log_file(const char * newname)
{
	PROFILE_FUNC();
	flush();
	if(m_fileStream.is_open()){
		UG_LOG("Logfile '" << m_logFileName << "' renamed to '" << newname << "'" << std::endl);
		rename(m_logFileName.c_str(), newname);
	}
	m_logFileName = newname;
	return true;
}

bool LogAssistant::
enable_terminal_output(bool bEnable)
{
	PROFILE_FUNC();
	m_terminalOutputEnabled = bEnable;
	update_ostream();
	return true;
}

void LogAssistant::
update_ostream()
{
	if(is_output_process()){
		if(m_terminalOutputEnabled){
			if(m_fileOutputEnabled){
			//	make sure that the logfile is enabled
				if(open_logfile())
					cout.rdbuf(m_splitBuf);
				else
					cout.rdbuf(m_logBuf);
			}
			else
				cout.rdbuf(m_logBuf);
		}
		else{
			if(m_fileOutputEnabled){
				if(open_logfile())
					cout.rdbuf(m_fileBuf);
				else
					cout.rdbuf(m_emptyBuf);
			}
			else{
				cout.rdbuf(m_emptyBuf);
			}
		}
	}
	else
		cout.rdbuf(m_emptyBuf);
}

bool LogAssistant::
open_logfile()
{
	PROFILE_FUNC();
	// open log file only when output process
	// otherwise slowdown on > 1024 cores.
	if(!is_output_process())
		return true;
	if(!m_fileStream.is_open()){
		m_fileStream.open(m_logFileName.c_str());
		if(!m_fileStream){
			m_fileOutputEnabled = false;
			UG_LOG("ERROR: LogAssistant::open_logfile failed: Couldn't open "
					<< m_logFileName << endl);
			return false;
		}
	}
	return true;
}


void LogAssistant::
flush_error_log()
{
	string str = m_errStream.str();
	if(!str.empty()){
		logger() << endl
				<< "********************************************************************************\n"
				<< "ERRORS OCCURED: " << endl << str << endl
				<< "********************************************************************************\n";
		m_errStream.clear();
	}
}




void LogAssistant::
set_output_process(int procRank)
{
	if(m_outputProc == procRank)
		return;
	
	int numProcs = 1;
	#ifdef UG_PARALLEL
		numProcs = pcl::NumProcs();
	#endif

	if(procRank < numProcs){
		m_outputProc = procRank;
		update_ostream();
	}
	else{
		string strProc;
		if(numProcs == 1)
			strProc = " process is ";
		else
			strProc = " processes are ";

		UG_LOG("WARNING: Won't change output process to " << procRank <<
				", since only " << numProcs << strProc << "available. "
				"Output process is left at " << m_outputProc << "." << endl);
	}
}

bool LogAssistant::
is_output_process()
{
	#ifdef UG_PARALLEL
		if((m_outputProc == pcl::ProcRank()) || (m_outputProc == -1))
			return true;
		else
			return false;
	#endif

	return true;
}

int LogAssistant::
get_process_rank()
{
	#ifdef UG_PARALLEL
		return pcl::ProcRank();
	#endif

	return 0;
}

int LogAssistant::get_debug_level_noninline(const char *debugID) const
{
	return GetDebugIDManager().get_debug_level(debugID);
}

bool LogAssistant::set_debug_level_noninline(const char* debugID, int level)
{
	return GetDebugIDManager().set_debug_level(debugID, level);
}



}//	end of namespace
