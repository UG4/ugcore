/*
 * log.cpp
 *
 *  Created on: 10.03.2010
 *      Author: andreasvogel, sebastianreiter
 */
#include <iostream>
#include "common/log.h"
using namespace std;

namespace ug{

LogAssistant::
LogAssistant() :
	m_terminalOutputEnabled(true),
	m_fileOutputEnabled(false)
{
}

LogAssistant::LogAssistant(const LogAssistant&)
{
}

void LogAssistant::init()
{
//	originally this was placed in the constructor.
//	However, an internal compiler error in MinGW (windows)
//	appears then. Having moved the stuff here into this init
//	method seems to have helped.
	set_debug_levels(-1);
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

bool LogAssistant::
enable_file_output(bool bEnable, const char* filename)
{
	if(bEnable){
		if(!m_fileStream.is_open())
		  {
			m_fileStream.open(filename);
			if(!m_fileStream)
			{
				m_fileOutputEnabled = false;
				return false;
			}
			m_logFileName = filename;
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
		if(!m_fileStream.is_open()){
			UG_LOG("Unable to rename logfile to '" << newname << "': no logfile open!");
			return false;
		} else {
			UG_LOG("Logfile '" << m_logFileName << "' renamed to '" << newname << "'" << std::endl);
			rename(m_logFileName, newname);
			m_logFileName = newname;

		}
	return true;
}

bool LogAssistant::
enable_terminal_output(bool bEnable)
{
	m_terminalOutputEnabled = bEnable;
	update_ostream();
	return true;
}

void LogAssistant::
update_ostream()
{
	if(m_terminalOutputEnabled){
		if(m_fileOutputEnabled)
			cout.rdbuf(m_splitBuf);
		else
			cout.rdbuf(m_logBuf);
	}
	else{
		if(m_fileOutputEnabled){
			cout.rdbuf(m_fileBuf);
		}
		else{
			cout.rdbuf(m_emptyBuf);
		}
	}
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
/*
std::ostream&
LogAssistant::
file_logger()
{
	static std::ofstream out("uglog.log");
	return out;
}
*/
}
